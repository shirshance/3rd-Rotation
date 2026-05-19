import os
import subprocess

output_folder = "results_TGC"
os.makedirs(output_folder, exist_ok=True)

codon = "TGC"
sites = ["emptyA", "Asite", "Psite"]
directions = ["P1", "M1", "M2"]

for site in sites:
    for direction in directions:
        output_name = os.path.join(output_folder, f"{codon}_{site}_{direction}")

        cmd = [
            "py",
            "mt_codon_translation_updated.py",
            "-f", "mouse_mt13.fasta",
            "-c", codon,
            "-s", site,
            "-di", direction,
            "-o", output_name,
            "-l", "OOF"
        ]

        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)

print("Done. 9 files were generated.")