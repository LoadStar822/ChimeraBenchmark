#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import defaultdict
from multitax import NcbiTx
import numpy as np
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


def read_result_file(filepath, standard_labels, software='default'):
    """
    Reads the result file and returns a dictionary mapping SequenceID to predicted TAXID.

    Parameters:
    - filepath: Path to the result file.
    - standard_labels: Dictionary of standard labels (sequence IDs).
    - software: The format of the result file (e.g., 'chimera', 'kraken2', 'ganon', 'ganon2', 'taxor').

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
                    # Kraken2 typically has a status character in the first column
                    # e.g., 'C' for classified, 'U' for unclassified
                    status = parts[0]
                    sequence_id = parts[1]
                    taxid = parts[2] if status == 'C' else 'unclassified'
                    predicted_labels[sequence_id] = taxid
                else:
                    # If the line doesn't have enough parts, mark as unclassified
                    sequence_id = parts[0]
                    predicted_labels[sequence_id] = 'unclassified'
            elif software in ['ganon', 'ganon2']:
                # Handle Ganon and Ganon2 result format
                # Example Ganon format: [sequence_id] \t [taxid] \t [other fields...]
                parts = line.split('\t')
                sequence_id = parts[0]
                if len(parts) >= 2:
                    taxid = parts[1]
                    if taxid == '-' or taxid == '':
                        predicted_labels[sequence_id] = 'unclassified'
                    else:
                        predicted_labels[sequence_id] = taxid
                else:
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
                    predicted_labels[sequence_id] = 'unclassified'
            else:
                # Handle unsupported software formats
                raise ValueError(f"Unsupported software format: {software}")
    return predicted_labels


def read_bracken_abundance(filepath, tax, rank):
    """
    Reads the Bracken output file and returns a dictionary mapping TAXID to abundance (fraction)
    at the specified taxonomic rank.

    Parameters:
    - filepath: Path to the Bracken result file.
    - tax: NcbiTx instance for taxonomy operations.
    - rank: Taxonomic rank to consider.

    Returns:
    - abundances: Dictionary with TAXID as keys and abundances (normalized fractions) as values.
    """
    abundances = defaultdict(float)
    total_fraction = 0.0
    with open(filepath, 'r', encoding='utf-8') as f:
        next(f)  # Skip header line
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 7:
                name = parts[0]
                taxid = parts[1]
                taxonomy_lvl = parts[2]
                fraction = float(parts[6])  # fraction_total_reads
                # Map the taxid to the specified rank
                try:
                    lineage = tax.lineage(taxid, ranks=[rank])
                    rank_taxid = lineage[0] if lineage else None
                    if rank_taxid:
                        abundances[rank_taxid] += fraction
                        total_fraction += fraction
                except KeyError:
                    continue
    # Normalize abundances
    if total_fraction > 0:
        for taxid in abundances:
            abundances[taxid] /= total_fraction
    return abundances


def compute_abundance(labels, tax, rank):
    """
    Calculate the abundance of a specified classification level based on the provided labels (standards or predictions).
    Return a dictionary with taxid as the key and abundance (proportion of sequences) as the value.
    """
    abundance = defaultdict(int)
    total_sequences = 0

    for seq_id, taxid in labels.items():
        if taxid == 'unclassified':
            continue
        try:
            lineage = tax.lineage(taxid, ranks=[rank])
            rank_taxid = lineage[0] if lineage else None
            if rank_taxid:
                abundance[rank_taxid] += 1
                total_sequences += 1
        except KeyError:
            continue

    for taxid in abundance:
        abundance[taxid] = abundance[taxid] / total_sequences if total_sequences > 0 else 0

    return abundance


def compute_l1_distance(abundance1, abundance2):
    """
    Compute the L1 distance between two abundance distributions.

    Parameters:
    - abundance1: Dictionary of abundances (taxid: abundance)
    - abundance2: Dictionary of abundances (taxid: abundance)

    Returns:
    - distance: L1 distance between the two distributions
    """
    all_taxids = set(abundance1.keys()).union(set(abundance2.keys()))
    distance = 0.0
    for taxid in all_taxids:
        freq1 = abundance1.get(taxid, 0)
        freq2 = abundance2.get(taxid, 0)
        distance += abs(freq1 - freq2)
    return distance


def evaluate_abundance(standard_file, result_files, software_list, output_file, dataset, database, tax=None):
    if tax is None:
        tax = NcbiTx()

    if len(result_files) != len(software_list):
        raise ValueError("The number of result files and software formats must match.")

    standard_labels = read_standard_file(standard_file)

    # Compute standard abundances once for all ranks
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    standard_abundances = {}
    for rank in ranks:
        abundance = compute_abundance(standard_labels, tax, rank)
        standard_abundances[rank] = abundance

    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Dataset Name', 'Database', 'Taxonomic Rank', 'Software', 'L1 Distance']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for filepath, software in zip(result_files, software_list):
            if software == 'bracken':
                # For Bracken, read abundances directly
                for rank in ranks:
                    predicted_abundance = read_bracken_abundance(filepath, tax, rank)
                    standard_abundance = standard_abundances[rank]
                    l1_distance = compute_l1_distance(standard_abundance, predicted_abundance)
                    writer.writerow({
                        'Dataset Name': dataset,
                        'Database': database,
                        'Taxonomic Rank': rank.capitalize(),
                        'Software': software,
                        'L1 Distance': f"{l1_distance:.4f}"
                    })
            else:
                predicted_labels = read_result_file(filepath, standard_labels, software=software)
                for rank in ranks:
                    predicted_abundance = compute_abundance(predicted_labels, tax, rank)
                    standard_abundance = standard_abundances[rank]
                    l1_distance = compute_l1_distance(standard_abundance, predicted_abundance)
                    writer.writerow({
                        'Dataset Name': dataset,
                        'Database': database,
                        'Taxonomic Rank': rank.capitalize(),
                        'Software': software,
                        'L1 Distance': f"{l1_distance:.4f}"
                    })

    print(f"Abundance analysis results have been saved to {output_file}.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Perform abundance analysis and compute L1 distance between standard and software results.')
    parser.add_argument('-s', '--standard', required=True, help='Path to the standard file')
    parser.add_argument('-r', '--results', required=True, nargs='+',
                        help='Paths to the result files from different software')
    parser.add_argument('-w', '--software', required=True, nargs='+',
                        help='Software formats corresponding to the result files (e.g., kraken2, ganon2, bracken)')
    parser.add_argument('-o', '--output', default='abundance_analysis.csv', help='Output CSV file name')
    parser.add_argument('-d', '--dataset', required=True, help='Name of the dataset')
    parser.add_argument('-db', '--database', required=True, help='Name of the database')
    args = parser.parse_args(argv)

    evaluate_abundance(
        standard_file=args.standard,
        result_files=args.results,
        software_list=args.software,
        output_file=args.output,
        dataset=args.dataset,
        database=args.database
    )


if __name__ == '__main__':
    main()
