from Bio import SeqIO
import matplotlib.pyplot as plt

target_gene = "ND4"

aa_codons = {
    "Serine": ["TCT", "TCC", "TCA", "TCG"],
    "Trp": ["TGG", "TGA"],      # mitochondrial Trp
    "Tyr": ["TAT", "TAC"],
    "Phe": ["TTT", "TTC"],
    "Cys": ["TGT", "TGC"]
}

# Build codon -> amino acid dictionary
codon_to_aa = {}

for aa, codon_list in aa_codons.items():
    for codon in codon_list:
        codon_to_aa[codon] = aa

# Set of all TRIT1-sensitive codons
trit1_codons = set(codon_to_aa.keys())

record = SeqIO.read("mouse_mtDNA.gb", "genbank")

for feature in record.features:

    if feature.type != "CDS":
        continue

    gene_name = feature.qualifiers["gene"][0].upper()

    if gene_name != target_gene:
        continue

    sequence = str(feature.extract(record.seq)).upper()

    codons = [
        sequence[i:i+3]
        for i in range(0, len(sequence) - 2, 3)
    ]

    protein_length = len(codons)

    plt.figure(figsize=(18, 3))

    # Main protein line
    plt.hlines(
        y=0,
        xmin=0,
        xmax=protein_length,
        color="black",
        linewidth=6
    )

    # Find consecutive TRIT1 codons
    for i in range(len(codons) - 1):

        codon1 = codons[i]
        codon2 = codons[i + 1]

        if codon1 in trit1_codons and codon2 in trit1_codons:

            aa1 = codon_to_aa[codon1]
            aa2 = codon_to_aa[codon2]

            # Blue line
            plt.vlines(
                x=i,
                ymin=-0.5,
                ymax=0.5,
                color="blue",
                linewidth=2
            )

            # Write amino acids above
            label = f"{aa1[:3]}-{aa2[:3]}"

            plt.text(
                i,
                0.7,
                label,
                rotation=90,
                fontsize=8,
                ha="center"
            )

            # Optional: print full info
            print(
                f"Position {i}: "
                f"{codon1}-{codon2} "
                f"({aa1}-{aa2})"
            )

    plt.title(f"{target_gene}: consecutive TRIT1-sensitive codon pairs")

    plt.xlabel("Codon position")
    plt.yticks([])

    plt.tight_layout()
    plt.show()

    break