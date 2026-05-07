from Bio import SeqIO
import matplotlib.pyplot as plt

# RNA codons converted to DNA codons
target_codons = [
    "TTT", "TTC",        # UUU, UUC = Phe
    "TAT", "TAC",        # UAU, UAC = Tyr
    "TGT", "TGC",        # UGU, UGC = Cys
    "TGG",               # UGG = Trp
    "TCT", "TCC", "TCA", "TCG"   # UCU, UCC, UCA, UCG = Ser(UCN)
]

record = SeqIO.read("mouse_mtDNA.gb", "genbank")

results = []

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0]
        sequence = str(feature.extract(record.seq)).upper()

        total_codons = len(sequence) // 3
        target_count = 0

        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3]

            if codon in target_codons:
                target_count += 1

        enrichment = (target_count / total_codons) * 100

        results.append((gene_name, enrichment, target_count, total_codons))

# Sort from most enriched to least enriched
results.sort(key=lambda x: x[1], reverse=True)

genes = [x[0] for x in results]
frequencies = [x[1] for x in results]

print("\nTRIT1-related codon enrichment:\n")

for gene, enrichment, count, total in results:
    print(f"{gene}: {count}/{total} codons ({enrichment:.2f}%)")

# Plot
plt.figure(figsize=(10, 5))

plt.bar(genes, frequencies)

plt.xlabel("Mitochondrial protein-coding gene")
plt.ylabel("TRIT1-related codons / total codons (%)")
plt.title("Enrichment of TRIT1-related codons across mtDNA-encoded proteins")

plt.xticks(rotation=45)

plt.tight_layout()
plt.show()