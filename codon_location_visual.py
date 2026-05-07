from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# UCU in RNA = TCT in DNA
target_codon = "TCT"
target_gene = "ND1"

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0].upper()

        if gene_name == target_gene:

            sequence = str(feature.extract(record.seq)).upper()
            gene_length = len(sequence)

            codon_positions = []

            # Check codon by codon
            for i in range(0, gene_length, 3):
                codon = sequence[i:i+3]

                if codon == target_codon:
                    codon_positions.append(i)

            # Create figure
            fig, ax = plt.subplots(figsize=(12, 2))

            # White gene background
            ax.add_patch(
                patches.Rectangle(
                    (0, 0), gene_length, 1,
                    facecolor="white",
                    edgecolor="black"
                )
            )

            # Blue codon regions
            for pos in codon_positions:
                ax.add_patch(
                    patches.Rectangle(
                        (pos, 0), 3, 1,
                        facecolor="blue",
                        edgecolor="blue"
                    )
                )

            ax.set_xlim(0, gene_length)
            ax.set_ylim(0, 1)
            ax.set_xlabel("Position in ND3 gene (bp)")
            ax.set_yticks([])
            ax.set_title(f"{target_gene} gene map: UCU codons")

            plt.show()

            print(f"{target_gene} length: {gene_length} bp")
            print(f"Number of UCU codons: {len(codon_positions)}")
            print("Positions:", codon_positions)

            break