import glob

output_file = "merged.fasta"

with open(output_file, "w") as outfile:
    for fasta in glob.glob("*.fasta"):   # takes all .fasta files in current folder
        print(f"Merging: {fasta}")

        with open(fasta) as infile:
            outfile.write(infile.read())

            # add newline between files
            outfile.write("\n")

print("Done! Created:", output_file)