from collections import Counter
import re

print("Paste OUTPUT 1, then press Enter twice:")
output1 = ""
while True:
    line = input()
    if line == "":
        break
    output1 += line + "\n"

print("Paste OUTPUT 2, then press Enter twice:")
output2 = ""
while True:
    line = input()
    if line == "":
        break
    output2 += line + "\n"

print("Paste OUTPUT 3, then press Enter twice:")
output3 = ""
while True:
    line = input()
    if line == "":
        break
    output3 += line + "\n"

all_outputs = output1 + output2 + output3

# Find codon pairs, e.g. TTC-TGA
codon_pairs = re.findall(r"\b[A-Z]{3}-[A-Z]{3}\b", all_outputs)

# Find amino acid pairs, e.g. Phe-Trp
aa_pairs = re.findall(r"\((.*?)\)", all_outputs)

codon_counter = Counter(codon_pairs)
aa_counter = Counter(aa_pairs)

print("\nMost common codon pairs:")
for pair, count in codon_counter.most_common():
    print(f"{pair}: {count}")

print("\nMost common amino acid pairs:")
for pair, count in aa_counter.most_common():
    print(f"{pair}: {count}")