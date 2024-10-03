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
    - software: The format of the result file (e.g., 'default', 'kraken2', 'ganon').

    Returns:
    - predicted_labels: Dictionary with SequenceID as keys and predicted TAXID as values.
    """
    predicted_labels = {seq_id: 'unclassified' for seq_id in standard_labels.keys()}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if software == 'chimera':
                # Handle chimera format
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
            elif software == 'ganon2':
                # Handle Ganon result format
                # Example Ganon format: [sequence_id] \t [taxid] \t [other fields...]
                parts = line.split('\t')
                if len(parts) >= 2:
                    sequence_id = parts[0]
                    taxid = parts[1]
                    predicted_labels[sequence_id] = taxid
            else:
                # Handle unsupported software formats
                raise ValueError(f"Unsupported software format: {software}")
    return predicted_labels


def compute_abundance(labels, tax, rank):
    """
    Calculate the abundance of a specified classification level based on the provided labels (standards or predictions).
    Return a dictionary with taxid as the key and abundance (number of sequences) as the value.
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

    # 获取所有taxid的集合
    all_taxids = set(abundance1.keys()).union(set(abundance2.keys()))
    distance = 0.0
    for taxid in all_taxids:
        freq1 = abundance1.get(taxid, 0)
        freq2 = abundance2.get(taxid, 0)
        distance += abs(freq1 - freq2)
    return distance


def main():
    parser = argparse.ArgumentParser(
        description='Perform abundance analysis and compute L1 distance between standard and software results.')
    parser.add_argument('-s', '--standard', required=True, help='Path to the standard file')
    parser.add_argument('-r', '--results', required=True, nargs='+',
                        help='Paths to the result files from different software')
    parser.add_argument('-w', '--software', required=True, nargs='+',
                        help='Software formats corresponding to the result files (e.g., default, kraken2, ganon2)')
    parser.add_argument('-o', '--output', default='abundance_analysis.csv', help='Output CSV file name')
    parser.add_argument('-d', '--dataset', required=True, help='Name of the dataset')
    parser.add_argument('-db', '--database', required=True, help='Name of the database')
    args = parser.parse_args()

    if len(args.results) != len(args.software):
        raise ValueError("The number of result files and software formats must match.")

    tax = NcbiTx()

    standard_labels = read_standard_file(args.standard)

    predicted_labels_list = []
    for filepath, software in zip(args.results, args.software):
        predicted_labels = read_result_file(filepath, standard_labels, software=software)
        predicted_labels_list.append(predicted_labels)

    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    with open(args.output, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Dataset Name', 'Database', 'Taxonomic Rank', 'Software', 'L1 Distance']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        standard_abundances = {}
        for rank in ranks:
            abundance = compute_abundance(standard_labels, tax, rank)
            standard_abundances[rank] = abundance

        for predicted_labels, software in zip(predicted_labels_list, args.software):
            for rank in ranks:
                predicted_abundance = compute_abundance(predicted_labels, tax, rank)
                standard_abundance = standard_abundances[rank]
                l1_distance = compute_l1_distance(standard_abundance, predicted_abundance)
                writer.writerow({
                    'Dataset Name': args.dataset,
                    'Database': args.database,
                    'Taxonomic Rank': rank.capitalize(),
                    'Software': software,
                    'L1 Distance': f"{l1_distance:.4f}"
                })

    print(f"Abundance analysis results have been saved to {args.output}.")


if __name__ == '__main__':
    main()
