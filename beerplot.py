import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("three_target_codons_in_a_row.csv")

if df.empty:
    print("No hits found.")
    exit()

# Count how many times each triplet appears
df["count"] = df.groupby("codons")["codons"].transform("count")

plt.figure(figsize=(12, 6))

for codon_triplet in df["codons"].unique():
    sub = df[df["codons"] == codon_triplet]

    plt.scatter(
        sub["protein_position_start"],
        sub["gene"],
        s=80,
        label=codon_triplet
    )

plt.xlabel("Protein codon position")
plt.ylabel("Mitochondrial gene")
plt.title("3 consecutive target-codon clusters in mouse mitochondrial proteins")

# Put legend outside the plot
plt.legend(
    title="Codon triplet",
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
    fontsize=8
)

plt.tight_layout()
plt.savefig("clean_protein_map.png", dpi=300)
plt.show()