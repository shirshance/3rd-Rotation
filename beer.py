from Bio import SeqIO
import csv

INPUT_FILE = "mouse_mtDNA.gb"
OUTPUT_FILE = "three_target_codons_in_a_row.csv"

TARGET_CODONS = {
    "TCT", "TCC", "TCA", "TCG",
    "TGG",
    "TAT", "TAC",
    "TTT", "TTC",
    "TGT", "TGC"
}

record = SeqIO.read(INPUT_FILE, "genbank")

results = []

for feature in record.features:
    if feature.type != "CDS":
        continue

    gene = feature.qualifiers.get("gene", ["unknown"])[0]
    product = feature.qualifiers.get("product", ["unknown"])[0]

    # Extract the CDS sequence correctly, including complement genes like ND6
    cds_seq = str(feature.extract(record.seq)).upper()

    # Split CDS into codons
    codons = [
        cds_seq[i:i+3]
        for i in range(0, len(cds_seq) - 2, 3)
    ]

    # Scan 3 codons in a row
    for i in range(len(codons) - 2):
        triplet = codons[i:i+3]

        if all(codon in TARGET_CODONS for codon in triplet):
            results.append({
                "gene": gene,
                "product": product,
                "codon_number_start": i + 1,
                "codons": "-".join(triplet),
                "protein_position_start": i + 1
            })

# Print results
for r in results:
    print(
        f"{r['gene']} | codon {r['codon_number_start']} | "
        f"{r['codons']} | {r['product']}"
    )

# Save to CSV
with open(OUTPUT_FILE, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "gene",
            "product",
            "codon_number_start",
            "codons",
            "protein_position_start"
        ]
    )
    writer.writeheader()
    writer.writerows(results)

print(f"\nDone. Found {len(results)} hits.")
print(f"Saved results to: {OUTPUT_FILE}")