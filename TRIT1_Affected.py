from Bio import SeqIO
import matplotlib.pyplot as plt

target_gene = "COX1"

# Amino acid codons in DNA version
aa_codons = {
    "Serine": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Trp": ["TGG"],
    "Tyr": ["TAT", "TAC"],
    "Phe": ["TTT", "TTC"],
    "Cys": ["TGT", "TGC"]
}

# Load GenBank file
record = SeqIO.read("mouse_mtDNA.gb", "genbank")

for feature in record.features:

    if feature.type == "CDS":

        gene_name = feature.qualifiers["gene"][0].upper()

        if gene_name == target_gene:

            sequence = str(feature.extract(record.seq)).upper()
            total_codons = len(sequence) // 3

            counts = {}

            for aa in aa_codons:
                counts[aa] = 0

            # Count codons
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]

                for aa, codon_list in aa_codons.items():
                    if codon in codon_list:
                        counts[aa] += 1

            # Calculate frequencies
            labels = list(counts.keys())
            frequencies = []

            for aa in labels:
                freq = (counts[aa] / total_codons) * 100
                frequencies.append(freq)

            # Print results
            print(f"\nAmino acid codon frequency in {target_gene}:\n")
            print(f"Total codons: {total_codons}\n")

            for aa in labels:
                print(f"{aa}: {counts[aa]} codons ({(counts[aa] / total_codons) * 100:.2f}%)")

            # Plot
            plt.figure(figsize=(8, 5))
            plt.bar(labels, frequencies)

            plt.xlabel("Amino acid")
            plt.ylabel("Frequency among all codons (%)")
            plt.title(f"{target_gene}: Ser/Trp/Tyr/Phe/Cys codon frequency")

            plt.tight_layout()
            plt.show()

            break