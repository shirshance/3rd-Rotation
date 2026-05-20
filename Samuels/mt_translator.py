# translate_dna.py

from mt_codon_table import CODON_TABLE

dna = input("Enter DNA sequence: ").upper().replace(" ", "")

protein = ""

for i in range(0, len(dna) - 2, 3):
    codon = dna[i:i+3]

    # Skip incomplete codons at the end
    if len(codon) != 3:
        continue

    aa = CODON_TABLE.get(codon, "X")  # X = unknown codon
    protein += aa

print("\nProtein:")
print(protein)