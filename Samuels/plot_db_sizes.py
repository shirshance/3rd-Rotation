#!/usr/bin/env python3
"""
Plot DB size comparisons from the CSV produced by compute_db_sizes.py

Generates two barplots (PNG):
- num_unique_sequences.png
- total_aa.png

Usage:
  python plot_db_sizes.py --csv db_sizes.csv
"""
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', required=False, default='db_sizes.csv')
    parser.add_argument('--outdir', default='plots')
    args = parser.parse_args()

    here = os.path.abspath(os.path.dirname(__file__))
    csv_path = os.path.join(here, args.csv)
    if not os.path.isfile(csv_path):
        raise SystemExit(f'CSV not found: {csv_path} — run compute_db_sizes.py first')

    df = pd.read_csv(csv_path)
    df = df.sort_values('num_unique_sequences', ascending=False)

    outdir = os.path.join(here, args.outdir)
    os.makedirs(outdir, exist_ok=True)

    sns.set(style='whitegrid')

    plt.figure(figsize=(10,6))
    ax = sns.barplot(x='db_name', y='num_unique_sequences', data=df, palette='viridis')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel('Number of unique sequences')
    plt.tight_layout()
    out1 = os.path.join(outdir, 'num_unique_sequences.png')
    plt.savefig(out1, dpi=150)
    plt.close()

    plt.figure(figsize=(10,6))
    df2 = df.sort_values('total_aa', ascending=False)
    ax = sns.barplot(x='db_name', y='total_aa', data=df2, palette='magma')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_ylabel('Total amino acids (unique sequences)')
    plt.tight_layout()
    out2 = os.path.join(outdir, 'total_aa.png')
    plt.savefig(out2, dpi=150)
    plt.close()

    print('Wrote plots to', outdir)


if __name__ == '__main__':
    main()
