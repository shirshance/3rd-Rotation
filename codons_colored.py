from Bio import SeqIO
import matplotlib.pyplot as plt

target_gene = "ND4"

# TRIT1-related codons + colors
aa_codons = {
    "Serine": {
        "codons": ["TCT", "TCC", "TCA", "TCG"],
        "color": "green"
    },
    "Trp": {
        "codons": ["TGG", "TGA"],   # mitochondrial Trp
        "color": "blue"
    },
    "Tyr": {
        "codons": ["TAT", "TAC"],
        "color": "orange"
    },
    "Phe": {
        "codons": ["TTT", "TTC"],
        "color": "purple"
    },
    "Cys": {
        "codons": ["TGT", "TGC"],
        "color": "pink"
    }
}

# Load mtDNA GenBank
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

for feature in record.features:

    if feature.type != "CDS":
        continue

    gene_name = feature.qualifiers["gene"][0].upper()

    if gene_name != target_gene:
        continue

    sequence = str(feature.extract(record.seq)).upper()

    # Split into codons
    codons = [
        sequence[i:i+3]
        for i in range(0, len(sequence), 3)
    ]

    protein_length = len(codons)

    # Create plot
    plt.figure(figsize=(16, 2))

    # Main protein line
    plt.hlines(
        y=0,
        xmin=0,
        xmax=protein_length,
        color="black",
        linewidth=6
    )

    # Plot codon positions
    for i, codon in enumerate(codons):

        for aa, info in aa_codons.items():

            if codon in info["codons"]:

                plt.vlines(
                    x=i,
                    ymin=-0.5,
                    ymax=0.5,
                    color=info["color"],
                    linewidth=2,
                    label=aa
                )

    # Clean duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()

    unique = dict(zip(labels, handles))

    plt.legend(
        unique.values(),
        unique.keys(),
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.title(f"{target_gene} TRIT1-sensitive codon map")

    plt.xlabel("Codon position")
    plt.yticks([])

    plt.tight_layout()
    plt.show()

    break