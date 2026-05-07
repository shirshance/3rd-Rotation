from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# UCU in RNA = TCT in DNA
target_pair = ["TCT", "TCT"]

record = SeqIO.read("mouse_mtDNA.gb", "genbank")

genes_data = []

for feature in record.features:
    if feature.type == "CDS":
        gene_name = feature.qualifiers["gene"][0]
        sequence = str(feature.extract(record.seq)).upper()
        gene_length = len(sequence)

        pair_positions = []

        # check codon by codon
        for i in range(0, gene_length - 3, 3):
            codon1 = sequence[i:i+3]
            codon2 = sequence[i+3:i+6]

            if [codon1, codon2] == target_pair:
                pair_positions.append(i)

        genes_data.append((gene_name, gene_length, pair_positions))

# Sort genes alphabetically
genes_data.sort(key=lambda x: x[0])

fig, ax = plt.subplots(figsize=(14, 8))

y = 0

for gene_name, gene_length, pair_positions in genes_data:
    # white gene bar
    ax.add_patch(
        patches.Rectangle(
            (0, y), gene_length, 0.6,
            facecolor="white",
            edgecolor="black"
        )
    )

    # blue marks for UCU-UCU pairs
    for pos in pair_positions:
        ax.add_patch(
            patches.Rectangle(
                (pos, y), 6, 0.6,
                facecolor="blue",
                edgecolor="blue"
            )
        )

    ax.text(-120, y + 0.3, gene_name, va="center", ha="right")
    y += 1

ax.set_xlim(-200, max(g[1] for g in genes_data) + 100)
ax.set_ylim(-0.5, y)
ax.set_xlabel("Position within gene (bp)")
ax.set_yticks([])
ax.set_title("UCU-UCU consecutive codons across mouse mtDNA protein-coding genes")

plt.tight_layout()
plt.show()

print("UCU-UCU pairs per gene:\n")
for gene_name, gene_length, pair_positions in genes_data:
    print(f"{gene_name}: {len(pair_positions)} pairs")