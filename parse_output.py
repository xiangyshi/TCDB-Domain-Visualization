import os
import re

# Directory containing the files to scan
RESCUED_DIR = "rescued"

# Function to check if any field in a line has more than 2 "|" characters
def has_too_many_pipes(line):
    fields = line.strip().split('\t')
    for field in fields:
        if field.count('|') > 2:
            return True
    return False

# Keep track of files with issues
files_with_issues = set()

# Scan all files in the rescued directory
for filename in os.listdir(RESCUED_DIR):
    filepath = os.path.join(RESCUED_DIR, filename)
    
    # Skip if not a file
    if not os.path.isfile(filepath):
        continue
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comments (lines starting with #)
            if line.startswith('#'):
                continue
            
            # Check if any field has more than 2 "|" characters
            if has_too_many_pipes(line):
                files_with_issues.add(filename)
                break  # No need to check more lines in this file

# Print the results
print(f"Files containing fields with more than 2 '|' characters:")
for filename in sorted(files_with_issues):
    print(f"- {filename}")
