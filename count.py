import os
for filename in os.listdir("holes"):
    if filename.endswith(".fasta"):
        with open("holes/" + filename, "r") as file:
            if len(file.readlines()) >= 20:
                print(filename)