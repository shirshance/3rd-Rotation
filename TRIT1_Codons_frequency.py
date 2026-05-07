from Bio import SeqIO

# Codons of interest
target_codons = [
    "TCT", "TCC", "TCA", "TCG",
    "TGG",
    "TAT", "TAC",
    "TTT", "TTC",
    "TGT", "TGC"
]

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

results = []

# Go over all features
for feature in record.features:

    # Only protein-coding genes
    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0]

        # Extract DNA sequence
        sequence = str(feature.extract(record.seq)).upper()

        total_bp = len(sequence)
        total_codons = total_bp // 3

        match_count = 0

        # Read codons
        for i in range(0, len(sequence), 3):

            codon = sequence[i:i+3]

            if codon in target_codons:
                match_count += 1

        # Frequency
        frequency = (match_count / total_codons) * 100

        results.append((gene_name, frequency, match_count, total_codons))

# Sort from highest to lowest frequency
results.sort(key=lambda x: x[1], reverse=True)

print("\nTRIT1-related codon frequencies:\n")

for gene, freq, matches, total in results:

    print(f"{gene}: {matches}/{total} codons ({freq:.2f}%)")