from Bio import SeqIO
import matplotlib.pyplot as plt

target_gene = "COX1"

# Standard mitochondrial amino acid codons (DNA version)
aa_codons = {
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Asn": ["AAT", "AAC"],
    "Asp": ["GAT", "GAC"],
    "Cys": ["TGT", "TGC"],
    "Gln": ["CAA", "CAG"],
    "Glu": ["GAA", "GAG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "His": ["CAT", "CAC"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Lys": ["AAA", "AAG"],
    "Met": ["ATG"],
    "Phe": ["TTT", "TTC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Trp": ["TGG"],
    "Tyr": ["TAT", "TAC"],
    "Val": ["GTT", "GTC", "GTA", "GTG"]
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

            # Frequencies
            labels = list(counts.keys())
            frequencies = []

            for aa in labels:

                freq = (counts[aa] / total_codons) * 100
                frequencies.append(freq)

            # Print results
            print(f"\nAmino acid frequencies in {target_gene}:\n")

            for aa in labels:
                print(f"{aa}: {counts[aa]} codons ({(counts[aa] / total_codons) * 100:.2f}%)")

            # Plot
            plt.figure(figsize=(12, 6))

            plt.bar(labels, frequencies)

            plt.xlabel("Amino acid")
            plt.ylabel("Frequency among all codons (%)")
            plt.title(f"{target_gene} amino acid composition")

            plt.xticks(rotation=45)

            plt.tight_layout()
            plt.show()

            break