#!/usr/bin/env python3

import os
import re
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import BiopythonWarning
import warnings
from datetime import datetime
from codon_table import CODON_TABLE


class NoNamePeptides:
    def __init__(self):
        self.noname_peptides = dict()
        self.named_peptides = set()

    def add(self, header, seq):
        if re.search(r'_NoName_', header):
            if seq in self.noname_peptides.keys():
                self.noname_peptides[seq].add(header)
            else:
                self.noname_peptides[seq] = {header}
        else:
            self.named_peptides.add(seq)

    def get(self):
        pept_array = []
        for seq in self.noname_peptides.keys():
            if seq not in self.named_peptides:
                for header in self.noname_peptides[seq]:
                    pept_array.append([header, seq])
        return pept_array


def validate_codons(codon):
    """
    Validate if the provided codon is a valid DNA codon (i.e., A, T, C, G).
    """
    if 'U' in codon:
        raise ValueError(f"Error: Codon '{codon}' contains RNA base 'U'. Please provide a DNA codon (A, T, C, G).")


def translate_shift(seq, description, cod, upstream, length, site, direction, stop_aa, table, trim):
    shifted = {}
    protein = ''
    count = 0

    for i in range(0, len(seq) - (len(seq) % 3), 3):
        codon = seq[i:i + 3]
        aa = table.get(codon, '_')  # Use default '_' if codon is not found in the table
        protein += aa

        # Debug: Track the sequence and amino acid translation
        print(f"Codon: {codon} -> Amino acid: {aa}")

        if aa == '_':  # Stop codon handling
            print(f"Stop codon encountered at position {i}. Ending protein translation.")
            break

        if codon == cod:  # Frameshift happens here
            count += 1
            print(f"Frameshift at position {i} with codon {codon} (Count: {count})")

            protein_in, protein_out = '', ''
            if site == 'emptyA':
                protein_in = protein[:-1]
                if direction == 'M1':
                    protein_out = str(seq[i - 1:].translate())
                elif direction == 'P1':
                    protein_out = str(seq[i + 1:].translate())
                elif direction == 'M2':
                    protein_out = str(seq[i - 2:].translate())

            elif site == 'Asite':
                protein_in = protein
                if direction == 'M1':
                    protein_out = str(seq[i + 3 - 1:].translate())
                elif direction == 'P1':
                    protein_out = str(seq[i + 3 + 1:].translate())
                elif direction == 'M2':
                    protein_out = str(seq[i + 3 - 2:].translate())

            elif site == 'Psite':
                next_aa = str(seq[i + 3:i + 6].translate())
                if next_aa == '*':  # The next codon can be a stop codon
                    continue
                protein_in = protein + next_aa
                if direction == 'M1':
                    protein_out = str(seq[i + 6 - 1:].translate())
                elif direction == 'P1':
                    protein_out = str(seq[i + 6 + 1:].translate())
                elif direction == 'M2':
                    protein_out = str(seq[i + 6 - 2:].translate())

            if len(protein_in) > int(upstream):
                protein_in = protein_in[-int(upstream):]

            protein_out = protein_out.split('*')[0]

            if stop_aa:
                protein_out = split_at_target_aa(protein_out, cod, aa, trim)

            # Debug: Show the protein_in and protein_out before writing
            print(f"Protein in: {protein_in}")
            print(f"Protein out: {protein_out}")

            identification_minus = f">{direction}_{site}_count{count}_{CODON_TABLE[cod]}{str(int(i / 3) + 1)}_" \
                                   f"{description.split(' ')[0]}_{protein_out[:14]}_{str(seq[i - 4:i + 6])}"

            shifted[identification_minus] = protein_in + protein_out

    return shifted


def split_at_target_aa(protein_out, cod, aa, trim):
    target_aa = CODON_TABLE[cod]
    protein_out = protein_out.split(target_aa)[0]
    if trim:
        protein_out = protein_out[:protein_out.find(target_aa)+1]
    return protein_out


def read_cds(fasta, cod, upstream, length, site, direction, stop_aa, gene_list, table, path_writing, trim):
    output_file = path_writing + '.fasta'

    entries = []  # collect (header, seq) tuples
    nn_peptide = NoNamePeptides()
    shifted_found = False

    for record in SeqIO.parse(fasta, "fasta"):
        seq = record.seq
        description = record.description
        gene = extract_gene_symbol(description)

        # Debugging: Print sequence and gene info
        print(f"Processing sequence: {description}")
        print(f"Sequence: {seq}")

        if gene_list:
            genes = pd.read_csv(gene_list)['Gene'].tolist()
            if gene not in genes:
                print(f"Skipping gene {gene} as it's not in the gene list.")
                continue  # Skip the gene if it's not in the list

        # Check if codon is present in the sequence
        if cod not in seq:
            print(f"Skipping {description} as codon {cod} is not in the sequence.")
            continue

        shifted = translate_shift(seq, description, cod, upstream, length, site, direction, stop_aa, table, trim)

        for k, v in shifted.items():
            entries.append((k, v))
            nn_peptide.add(k, v)
            shifted_found = True  # Mark that we've processed a peptide

    # Add any NoName peptides that are not duplicates
    for k, v in nn_peptide.get():
        entries.append((k, v))

    # Filter out peptides shorter than 8 amino acids
    min_len = 8
    retained = [(h, s) for (h, s) in entries if len(s) >= min_len]
    removed_count = len(entries) - len(retained)

    # Write filtered results to file
    with open(output_file, "w+") as out_file:
        for h, s in retained:
            out_file.write(h + "\n")
            out_file.write(s + "\n")

    # Reporting
    if removed_count > 0:
        print(f"Removed {removed_count} peptides shorter than {min_len} AA.")

    print(f"Wrote {len(retained)} peptides to {output_file}")

    if not shifted_found:
        print("Warning: No valid peptides were processed. Please check your input data.")

    if len(retained) == 0:
        print("Warning: No peptides met the minimum length threshold. Output FASTA is empty.")

    return


def extract_gene_symbol(description):
    if re.search(r'gene_symbol:', description):
        return str(description.split('gene_symbol:')[1].split(' ')[0])
    else:
        return 'NoName'


def create_output_filename(args):
    now = datetime.now()
    date = now.strftime("%m-%d-%Y")

    # Format string without references to `args.data`
    return '{}_{}_{}_{}_{}{}{}_{}'.format(
        args.output,                    # args.output is the filename prefix
        str(args.upstream),              # Upstream length
        str(args.codon),                 # Codon (target codon for frameshift)
        args.length,                     # Translation type (chimeras or OOF)
        ('_Stop' if args.stop_aa else '_NoStop'),  # Stop amino acid option
        CODON_TABLE.get(args.codon, 'UNKNOWN'),  # Get codon translation or 'UNKNOWN'
        ('Trimmed_' if args.trim else ''),  # Trimmed option
        date  # Date
    )


def main():
    parser = argparse.ArgumentParser(description="Frameshift Peptide Generator")
    parser.add_argument('-f', '--fasta', required=True, help='Path to FASTA file with CDS sequences')
    parser.add_argument('-c', '--codon', required=True, help='Target codon for frameshift')
    parser.add_argument('-s', '--site', choices=['emptyA', 'Asite', 'Psite'], help='Site type for frameshift')
    parser.add_argument('-di', '--direction', choices=['P1', 'M1', 'M2'], help='Direction of frameshift')
    parser.add_argument('-u', '--upstream', default=13, type=int, help='Upstream length')
    parser.add_argument('-l', '--length', choices=['chimeras', 'OOF'], default='chimeras', help='Translation type')
    parser.add_argument('-st', '--stop_aa', action='store_true', help='Stop at a target amino acid')
    parser.add_argument('-tr', '--trim', action='store_true', help='Trim by target amino acid')
    parser.add_argument('-g', '--gene_list', help="Path to gene list CSV")
    parser.add_argument('-o', '--output', required=True, help='Output file prefix')
    args = parser.parse_args()

    output_path = create_output_filename(args)

    # Now call the main translation code (no file paths here)
    read_cds(args.fasta, args.codon, args.upstream, args.length, 
             args.site, args.direction, args.stop_aa, args.gene_list, CODON_TABLE, output_path, args.trim)


if __name__ == '__main__':
    main()
