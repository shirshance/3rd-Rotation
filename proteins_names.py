from Bio import SeqIO

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

genes = []

# Go over features
for feature in record.features:

    # Only CDS features
    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0]

        # Extract sequence
        sequence = feature.extract(record.seq)

        # Sequence length
        length = len(sequence)

        genes.append((gene_name, length))

# Sort alphabetically
genes.sort()

print("Protein-coding mitochondrial genes:\n")

for gene, length in genes:
    print(f"{gene} - {length} bp")