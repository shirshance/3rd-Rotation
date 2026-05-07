from Bio import SeqIO

# Load the GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

# Ask user for gene name
gene_name = input("Enter mitochondrial gene name: ").upper()

found = False

# Go over all annotated features
for feature in record.features:

    # Check only genes/CDS
    if feature.type == "CDS":

        # Some features contain gene names
        if "gene" in feature.qualifiers:

            gene = feature.qualifiers["gene"][0].upper()

            # If user gene matches
            if gene == gene_name:

                # Extract DNA sequence
                sequence = feature.extract(record.seq)

                print(f"\nGene: {gene}")
                print(f"Sequence:\n{sequence}")

                found = True
                break

if not found:
    print("Gene not found.")
