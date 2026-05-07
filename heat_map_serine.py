from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# DNA codons -> RNA codon names
serine_codons = {
    "TCT": "UCU",
    "TCC": "UCC",
    "TCA": "UCA",
    "TCG": "UCG",
    "AGC": "AGC",
    "AGT": "AGU"
}

codon_labels = ["UCU", "UCC", "UCA", "UCG", "AGC", "AGU"]

record = SeqIO.read("mouse_mtDNA.gb", "genbank")

data = []
gene_names = []

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0]
        sequence = str(feature.extract(record.seq)).upper()

        counts = {codon: 0 for codon in codon_labels}

        # count serine codons
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]

            if codon in serine_codons:
                rna_codon = serine_codons[codon]
                counts[rna_codon] += 1

        total_serines = sum(counts.values())

        frequencies = []

        for codon in codon_labels:
            if total_serines > 0:
                freq = (counts[codon] / total_serines) * 100
            else:
                freq = 0

            frequencies.append(freq)

        gene_names.append(gene_name)
        data.append(frequencies)

data = np.array(data)

# Plot heatmap
plt.figure(figsize=(9, 7))

plt.imshow(data, aspect="auto")

plt.xticks(range(len(codon_labels)), codon_labels)
plt.yticks(range(len(gene_names)), gene_names)

plt.xlabel("Serine codon")
plt.ylabel("Mitochondrial gene")
plt.title("Serine codon usage across mouse mtDNA proteins")

plt.colorbar(label="Frequency among serines (%)")

# Add numbers inside cells
for i in range(len(gene_names)):
    for j in range(len(codon_labels)):
        plt.text(j, i, f"{data[i, j]:.1f}", ha="center", va="center")

plt.tight_layout()
plt.show()