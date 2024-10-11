#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from multitax import NcbiTx
import csv


def read_standard_file(filepath):
    """
    Reads the standard file and returns a dictionary mapping SequenceID to TAXID.

    Parameters:
    - filepath: Path to the standard file.

    Returns:
    - standard_labels: Dictionary with SequenceID as keys and TAXID as values.
    """
    standard_labels = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Skip header lines that start with '@'
    data_lines = [line for line in lines if not line.startswith('@')]

    for line in data_lines:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 3:
            sequence_id = parts[0]
            taxid = parts[2]
            standard_labels[sequence_id] = taxid
    return standard_labels


def read_result_file(filepath, standard_labels, software):
    """
    Reads the result file and returns a dictionary mapping SequenceID to predicted TAXID.

    Parameters:
    - filepath: Path to the result file.
    - standard_labels: Dictionary of standard labels (sequence IDs).
    - software: The format of the result file (e.g., 'chimera', 'kraken2', 'ganon', 'taxor').

    Returns:
    - predicted_labels: Dictionary with SequenceID as keys and predicted TAXID as values.
    """
    predicted_labels = {seq_id: 'unclassified' for seq_id in standard_labels.keys()}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if software == 'chimera':
                # Handle Chimera format
                parts = line.split('\t')
                sequence_id = parts[0]
                if len(parts) >= 2:
                    if parts[1].lower() == 'unclassified':
                        predicted_labels[sequence_id] = 'unclassified'
                    else:
                        # Extract the first predicted TAXID
                        taxid_score = parts[1]
                        taxid = taxid_score.split(':')[0]
                        predicted_labels[sequence_id] = taxid
                else:
                    predicted_labels[sequence_id] = 'unclassified'
            elif software == 'kraken2':
                # Handle Kraken2 result format
                # Example Kraken2 format: [status] \t [sequence_id] \t [taxonomy]
                parts = line.split('\t')
                if len(parts) >= 3:
                    # Kraken2 uses 'C' for classified, 'U' for unclassified
                    status = parts[0]
                    sequence_id = parts[1]
                    taxid = parts[2] if status == 'C' else 'unclassified'
                    predicted_labels[sequence_id] = taxid
                else:
                    # If the line doesn't have enough parts, mark as unclassified
                    sequence_id = parts[0]
                    predicted_labels[sequence_id] = 'unclassified'
            elif software in ['ganon', 'ganon2']:
                # Handle Ganon format (same for ganon and ganon2)
                # Example Ganon format: [sequence_id] \t [taxid] \t [other fields...]
                parts = line.split('\t')
                if len(parts) >= 2:
                    sequence_id = parts[0]
                    taxid = parts[1]
                    if taxid == '-':
                        predicted_labels[sequence_id] = 'unclassified'
                    else:
                        predicted_labels[sequence_id] = taxid
                else:
                    # If the line doesn't have enough parts, mark as unclassified
                    sequence_id = parts[0]
                    predicted_labels[sequence_id] = 'unclassified'
            elif software == 'taxor':
                # Handle Taxor format
                # Header: #QUERY_NAME	ACCESSION	REFERENCE_NAME	TAXID	REF_LEN	QUERY_LEN	QHASH_COUNT	QHASH_MATCH	TAX_STR	TAX_ID_STR
                parts = line.split('\t')
                sequence_id = parts[0]
                if len(parts) >= 4:
                    taxid = parts[3]
                    if taxid == '-' or taxid == '':
                        predicted_labels[sequence_id] = 'unclassified'
                    else:
                        predicted_labels[sequence_id] = taxid
                else:
                    # If the line doesn't have enough parts, mark as unclassified
                    predicted_labels[sequence_id] = 'unclassified'
            else:
                # Handle unsupported software formats
                raise ValueError(f"Unsupported software format: {software}")
    return predicted_labels


def compute_metrics_at_rank(standard_labels, predicted_labels, tax, rank):
    """
    Computes Accuracy, Precision, Recall, and F1 Score.

    Parameters:
    - y_true: List of true TAXIDs.
    - y_pred: List of predicted TAXIDs.

    Returns:
    - accuracy: Accuracy score.
    - precision: Precision score.
    - recall: Recall score.
    - f1: F1 score.
    """
    TP = 0
    FP = 0
    FN = 0
    total_samples = len(standard_labels)

    for seq_id in standard_labels:
        true_taxid = standard_labels[seq_id]
        pred_taxid = predicted_labels.get(seq_id, 'unclassified')

        if pred_taxid == 'unclassified':
            FN += 1
            continue

        try:
            true_lineage = tax.lineage(true_taxid, ranks=[rank])
            pred_lineage = tax.lineage(pred_taxid, ranks=[rank])
        except KeyError:
            FN += 1
            continue

        true_rank_taxid = true_lineage[0] if true_lineage else None
        pred_rank_taxid = pred_lineage[0] if pred_lineage else None

        if true_rank_taxid is None or pred_rank_taxid is None:
            FN += 1
        elif true_rank_taxid == pred_rank_taxid:
            TP += 1
        else:
            FP += 1

    accuracy = TP / total_samples if total_samples > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    if (precision + recall) > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = 0

    return TP, FP, FN, accuracy, precision, recall, f1


def evaluate_classification(standard_file, result_file, software, output_file, dataset, database, tax=None):
    if tax is None:
        tax = NcbiTx()
    # Read standard labels
    standard_labels = read_standard_file(standard_file)

    # Read predicted labels
    predicted_labels = read_result_file(result_file, standard_labels, software=software)

    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Dataset Name', 'Database', 'Taxonomic Rank', 'Total Samples', 'True Positives (TP)',
                      'False Positives (FP)', 'False Negatives (FN)', 'Accuracy', 'Precision', 'Recall', 'F1 Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for rank in ranks:
            TP, FP, FN, accuracy, precision, recall, f1 = compute_metrics_at_rank(standard_labels, predicted_labels,
                                                                                  tax, rank)
            writer.writerow({
                'Dataset Name': dataset,
                'Database': database,
                'Taxonomic Rank': rank.capitalize(),
                'Total Samples': len(standard_labels),
                'True Positives (TP)': TP,
                'False Positives (FP)': FP,
                'False Negatives (FN)': FN,
                'Accuracy': f"{accuracy:.4f}",
                'Precision': f"{precision:.4f}",
                'Recall': f"{recall:.4f}",
                'F1 Score': f"{f1:.4f}"
            })

    print(f"Evaluation metrics at different taxonomic levels have been saved to {output_file}.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Compare classification results with a standard and calculate evaluation metrics.')
    parser.add_argument('-s', '--standard', required=True, help='Path to the standard file')
    parser.add_argument('-r', '--result', required=True, help='Path to the result file')
    parser.add_argument('-w', '--software', help='Software format of the result file (e.g., kraken2, ganon2)')
    parser.add_argument('-o', '--output', default='classification_results.csv', help='Output CSV file name')
    parser.add_argument('-d', '--dataset', required=True, help='Name of the dataset')
    parser.add_argument('-db', '--database', required=True, help='Name of the database')
    args = parser.parse_args(argv)

    evaluate_classification(
        standard_file=args.standard,
        result_file=args.result,
        software=args.software,
        output_file=args.output,
        dataset=args.dataset,
        database=args.database
    )

if __name__ == '__main__':
    main()
