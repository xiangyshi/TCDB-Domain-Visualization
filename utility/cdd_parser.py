import subprocess
import platform
import pandas as pd
from config import *

def clean_cdd(in_file, process_fam_ids) -> pd.DataFrame:
    """
    Cleans CDD output data by removing comments and extracting relevant fields.

    Executes system-specific cleaning commands and converts the output into
    a structured DataFrame format.

    Args:
        in_file (str): Path to input CDD file
        process_fam_ids (list): List of TCIDs to process

    Returns:
        pd.DataFrame: Cleaned data with extracted fields as columns

    Raises:
        ValueError: If Fields header is not found in input file
        Exception: If cleaning command fails
    """

    # Adjust command pipeline based on machine architecture
    arch = platform.system()
    CLEAN_COMMAND = UNIX_CLEAN_COMMAND
    if arch == 'Windows':
        CLEAN_COMMAND = WIN_CLEAN_COMMAND
    
    # Append input_file path
    CLEAN_COMMAND += [in_file]

    try:
        # Obtain column headers
        header = None
        with open(in_file, 'r') as input_file:
            for line in input_file:
                if line.startswith("# Fields:"):
                    header = line.strip().replace("# Fields: ", "")  # Extract and clean the header line
                    break  # Exit after finding the header
        if not header:
            raise ValueError("Fields header not found in the input file.")
        
        # Run the CLEAN command
        clean_process = subprocess.run(CLEAN_COMMAND, capture_output=True, text=True)
        if clean_process.returncode != 0:
            raise Exception(f"Error in CLEAN command: {clean_process.stderr}")

        # Convert output to DataFrame
        from io import StringIO
        cleaned_data = StringIO(clean_process.stdout)
        df = pd.read_csv(cleaned_data, sep='\t', header=None)

        # Set column names
        df.columns = header.split(", ")

        # Extract family ID from query accession
        df['family'] = df['query acc.'].apply(lambda x: '.'.join(x.split(".")[:3]))
        
        # Filter by process_fam_ids if non-empty
        if process_fam_ids:
            df = df[df['family'].isin(process_fam_ids)]

        return df
    
    except Exception as e:
        # return empty df when error
        print("An error occurred during cleaning input/projection file:", e)
        return pd.DataFrame() 