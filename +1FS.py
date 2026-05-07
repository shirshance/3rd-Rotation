
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

gb_file = "mouse_mtDNA.gb"   # change if needed
record = SeqIO.read(gb_file, "genbank")

mt_genes = [
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8"
]

# Mitochondrial stop codons in DNA
# TGA is NOT stop in mitochondria; it codes for Trp
stop_codons = {"TAA", "TAG"}

results = []

for feature in record.features:
    if feature.type != "CDS":
        continue

    gene = feature.qualifiers.get("gene", [""])[0].upper()
    gene = gene.replace("MT-", "")

    if gene not in mt_genes:
        continue

    dna = str(feature.extract(record.seq)).upper()

    for shift_name, shifted_seq in {
        "+1": dna[1:],
        "-1": dna[2:]
    }.items():

        codons = [
            shifted_seq[i:i+3]
            for i in range(0, len(shifted_seq) - 2, 3)
        ]

        stop_positions = [
            i + 1
            for i, codon in enumerate(codons)
            if codon in stop_codons
        ]

        first_stop = stop_positions[0] if stop_positions else None

        results.append({
            "gene": gene,
            "complex": "Complex I" if gene.startswith("ND") else "Other",
            "shift": shift_name,
            "gene_length_codons": len(dna) // 3,
            "number_of_stop_codons": len(stop_positions),
            "first_stop_codon_position": first_stop,
            "stop_positions": stop_positions
        })

df = pd.DataFrame(results)

print(df[[
    "gene",
    "complex",
    "shift",
    "gene_length_codons",
    "number_of_stop_codons",
    "first_stop_codon_position"
]].sort_values(["shift", "first_stop_codon_position"]))

df.to_csv("frameshift_stop_analysis.csv", index=False)
print("\nSaved results to frameshift_stop_analysis.csv")




# Plot first stop codon position
plot_df = df.copy()

# If no stop was found, put it at the end of the gene
plot_df["first_stop_for_plot"] = plot_df["first_stop_codon_position"].fillna(
    plot_df["gene_length_codons"]
)

for shift in ["+1", "-1"]:
    sub = plot_df[plot_df["shift"] == shift].sort_values("first_stop_for_plot")

    plt.figure(figsize=(10, 5))
    plt.bar(sub["gene"], sub["first_stop_for_plot"])
    plt.ylabel("First stop codon position after frameshift")
    plt.xlabel("mtDNA gene")
    plt.title(f"Frameshift sensitivity: {shift} frame")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()