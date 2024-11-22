import os
import subprocess

# Directory containing your input files
input_dir = "holes"
# Directory where you want to save the output files
output_dir = "aligns"
# Scoring
score_dir = "scores"
print("Align in progress...")
# Loop through each .txt file in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):  # Ensure the file is a .fasta file
        input_file = os.path.join(input_dir, filename)
        output_file = os.path.join(output_dir, f"{os.path.splitext(filename)[0][:-6]}_align.aln")
        
        # Prepare the ClustalW command
        clustalw_command = [
            "clustalw2", 
            f"-INFILE={input_file}", 
            f"-OUTFILE={output_file}", 
            "-OUTPUT=CLUSTAL"
        ]
        
        # Run the ClustalW command using subprocess
        try:
            result = subprocess.run(clustalw_command, stdout=subprocess.PIPE, text=True)
            with open(os.path.join(score_dir, os.path.splitext(filename)[0][:-6] + '.log'), "w") as log_file:
                log_file.write(result.stdout)
            print(f"Alignment for {filename} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {filename}: {e}")
