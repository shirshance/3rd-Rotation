#!/usr/bin/env python3
"""Compute unique-sequence counts and total AA using an on-disk SQLite table.

Usage:
  python compute_sizes_by_sqlite.py path/to/file1.fasta [path/to/file2.fasta ...]
"""
import sys
import os
import sqlite3
from Bio import SeqIO


def compute_for_file(fpath):
    tmpdb = fpath + '.tmp.sqlite'
    if os.path.exists(tmpdb):
        os.remove(tmpdb)
    conn = sqlite3.connect(tmpdb)
    cur = conn.cursor()
    cur.execute('CREATE TABLE seqs(seq TEXT PRIMARY KEY)')
    conn.commit()
    inserted = 0
    for rec in SeqIO.parse(fpath, 'fasta'):
        s = str(rec.seq).strip()
        if not s:
            continue
        try:
            cur.execute('INSERT INTO seqs(seq) VALUES (?)', (s,))
            inserted += 1
        except sqlite3.IntegrityError:
            pass
        if inserted % 100000 == 0:
            conn.commit()
    conn.commit()
    cur.execute('SELECT COUNT(*) FROM seqs')
    num_unique = cur.fetchone()[0]
    # compute total aa from unique sequences
    total_aa = 0
    cur2 = conn.cursor()
    cur2.execute('SELECT seq FROM seqs')
    for (seq,) in cur2:
        total_aa += len(seq)
    conn.close()
    os.remove(tmpdb)
    return num_unique, total_aa


def main():
    if len(sys.argv) < 2:
        print('Usage: compute_sizes_by_sqlite.py file1.fasta [file2.fasta ...]')
        sys.exit(1)
    for f in sys.argv[1:]:
        if not os.path.isfile(f):
            print('Not found:', f)
            continue
        print('Processing', f)
        num_unique, total_aa = compute_for_file(f)
        avg = total_aa / num_unique if num_unique else 0
        print(os.path.basename(f), 'unique=', num_unique, 'total_aa=', total_aa, 'avg=', round(avg,2))


if __name__ == '__main__':
    main()
