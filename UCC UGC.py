from Bio import SeqIO

# UCC and UGC
# DNA equivalents:
# UCC -> TCC
# UGC -> TGC

target_codons = ["TCC", "TGC"]

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

results = []

# Go over all protein-coding genes
for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0]

        # Extract DNA sequence
        sequence = str(feature.extract(record.seq)).upper()

        total_codons = len(sequence) // 3

        match_count = 0

        # Read codons
        for i in range(0, len(sequence), 3):

            codon = sequence[i:i+3]

            if codon in target_codons:
                match_count += 1

        # Frequency
        frequency = (match_count / total_codons) * 100

        results.append((gene_name, frequency, match_count, total_codons))

# Sort by frequency
results.sort(key=lambda x: x[1], reverse=True)

print("\nUCC + UGC codon prevalence in mt genes:\n")

for rank, (gene, freq, matches, total) in enumerate(results, start=1):

    print(f"{rank}. {gene}: {matches}/{total} codons ({freq:.2f}%)")