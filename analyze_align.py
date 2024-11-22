import os
import shutil

def count(s):
    cnt_aster = 0
    for c in s:
        if c == "*":
            cnt_aster += 1
    return cnt_aster

# Directories
input_dir_1 = "holes"
input_dir_2 = "aligns"
important_dir = "important"

# Create the "important" directory if it doesn't exist
os.makedirs(important_dir, exist_ok=True)

lengths = {}
print("Analyzing...")

# Loop through each .fasta file in the input directory
for filename in os.listdir(input_dir_1):
    if filename.endswith(".fasta"):
        with open(os.path.join(input_dir_1, filename), "r") as file:
            lines = file.readlines()
            lengths[filename[:-12]] = (max([len(line) for line in lines]) - 1, len(lines) / 2)

# Loop through each .aln file in the input directory
for filename in os.listdir(input_dir_2):
    if filename.endswith(".aln"):
        with open(os.path.join(input_dir_2, filename), "r") as file:
            num_asters = count(file.read())
            percentage = num_asters / lengths[filename[:-10]][0]
            if percentage > 0.2 and lengths[filename[:-10]][1] > 2:
                # Copy the file to the "important" directory
                shutil.copy(os.path.join(input_dir_2, filename), important_dir)
                print(f"Copied: {filename} (Percentage: {percentage:.2f}, Length: {int(lengths[filename[:-10]][1])})")
