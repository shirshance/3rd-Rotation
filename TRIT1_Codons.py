from Bio import SeqIO

# Codons of interest
target_codons = [
    "TCT", "TCC", "TCA", "TCG",   # UCU UCC UCA UCG
    "TGG", "TGA",                     # UGG
    "TAT", "TAC",                 # UAU UAC
    "TTT", "TTC",                 # UUU UUC
    "TGT", "TGC"                  # UGU UGC
]

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

# Ask user for gene name
gene_input = input("Enter gene name: ").upper()

found = False

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0].upper()

        if gene_name == gene_input:

            # Extract DNA sequence
            sequence = str(feature.extract(record.seq)).upper()

            total_bp = len(sequence)
            total_codons = total_bp // 3

            match_count = 0

            # Read sequence codon by codon
            for i in range(0, len(sequence), 3):

                codon = sequence[i:i+3]

                if codon in target_codons:
                    match_count += 1

            # Frequency calculation
            frequency = (match_count / total_bp) * 100

            print(f"\nGene: {gene_name}")
            print(f"Total bp: {total_bp}")
            print(f"Total codons: {total_codons}")
            print(f"Matching codons: {match_count}")
            print(f"Frequency: {frequency:.2f}%")

            found = True
            break

if not found:
    print("Gene not found.")