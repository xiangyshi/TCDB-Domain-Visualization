import os

def read_asters(s):
    cnt = 0
    for c in s:
        if c == "*":
            cnt += 1
    return cnt

# Directory containing your input files
input_dir_1 = "holes"
input_dir_2 = "aligns"
lengths = {}
print("Analyzing...")
# Loop through each .txt file in the input directory
for filename in os.listdir(input_dir_1):
    if filename.endswith(".fasta"):
        with open(os.path.join(input_dir_1, filename), "r") as file:
            lines = file.readlines()
            lengths[filename[:-12]] = max([len(line) for line in lines]) - 1

for filename in os.listdir(input_dir_2):
    if filename.endswith(".aln"):
        with open(os.path.join(input_dir_2, filename), "r") as file:
            num_asters = read_asters(file.read())
            percentage = num_asters / lengths[filename[:-10]]
            if percentage > 0.1:
                print(filename, percentage)