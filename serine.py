from Bio import SeqIO
import matplotlib.pyplot as plt

target_gene = "ND3"

# Serine codons (DNA version)
serine_codons = {
    "TCT": "UCU",
    "TCC": "UCC",
    "TCA": "UCA",
    "TCG": "UCG",
    "AGC": "AGC",
    "AGT": "AGU"
}

# Initialize counts
counts = {
    "UCU": 0,
    "UCC": 0,
    "UCA": 0,
    "UCG": 0,
    "AGC": 0,
    "AGU": 0
}

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0].upper()

        if gene_name == target_gene:

            sequence = str(feature.extract(record.seq)).upper()

            # Go codon by codon
            for i in range(0, len(sequence), 3):

                codon = sequence[i:i+3]

                if codon in serine_codons:

                    rna_name = serine_codons[codon]
                    counts[rna_name] += 1

            break

# Total serines
total_serines = sum(counts.values())

# Frequencies
frequencies = []

labels = ["UCU", "UCC", "UCA", "UCG", "AGC", "AGU"]

for codon in labels:

    if total_serines > 0:
        freq = (counts[codon] / total_serines) * 100
    else:
        freq = 0

    frequencies.append(freq)

# Print results
print(f"\nSerine codon usage in {target_gene}:\n")

for codon in labels:
    print(f"{codon}: {counts[codon]}")

# Plot
plt.figure(figsize=(8,5))

plt.bar(labels, frequencies)

plt.xlabel("Serine codons")
plt.ylabel("Frequency among serines (%)")
plt.title(f"Serine codon usage in {target_gene}")

plt.show()