# TCDB Domain Visualization

A comprehensive Python tool for analyzing protein domain architectures from TCDB (Transport Classification Database) proteins using rpsblast and CDD (Conserved Domain Database) files.

## Overview

This project processes protein domains from TCDB-specific proteins, analyzing domain data from CDD files or rescue files to generate detailed domain architecture analysis and visualizations. It supports various visualization types including domain architecture plots, characteristic analysis, hole detection, and summary statistics.

## Features

- **Domain Architecture Analysis**: Extract and analyze protein domain arrangements
- **Multiple Input Formats**: Support for CDD files and rescue directories
- **Flexible Visualization**: Generate various plot types (architecture, characteristics, holes, summaries)
- **Data Export**: Export processed data to CSV format for further analysis
- **Targeted Processing**: Process specific TCID families or all available families
- **Overlapping Domain Handling**: Configurable merging of overlapping domain hits

## Installation

### Prerequisites

- Python 3.8 or higher
- All dependencies listed in `requirements.txt`

### Setup

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd TCDB-Domain-Visualization
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Ensure you have the necessary input data:
   - CDD files (e.g., `singles.cdd`)
   - Or rescue files in a directory structure
   - TCDB sequence files in `sequences/` directory

## Usage

The main CLI tool is `domain_extract.py`. Here's how to use it:

### Basic Command Structure

```bash
python domain_extract.py [INPUT_OPTIONS] [PROCESSING_OPTIONS] [OUTPUT_OPTIONS]
```

### Input Options (Required - Choose One)

#### Option 1: CDD File Input
```bash
python domain_extract.py -c path/to/your_file.cdd
```

#### Option 2: Rescue Directory Input
```bash
python domain_extract.py -r path/to/rescue_directory/
```

### Processing Options
#### Target Specific TCIDs
```bash
# Process specific families (comma-separated)
python domain_extract.py -c singles.cdd -t "1.A.12,1.C.105,2.A.1"

# Process families from a file (one TCID per line)
python domain_extract.py -c singles.cdd -t tcid_list.txt
```

### Output Options

#### Set Output Directory
```bash
python domain_extract.py -c singles.cdd -o results/
```

#### Generate Data Export
```bash
python domain_extract.py -c singles.cdd -d protein_data.csv
```

#### Generate Plots
```bash
# Generate all available plots
python domain_extract.py -c singles.cdd -p all

# Generate specific plot types for CDD input
python domain_extract.py -c singles.cdd -p "general,char,arch"

# Generate specific plot types for rescue input
python domain_extract.py -r rescue_folder/ -p "char_rescue"

# Available plot types:
# For CDD: general, char, arch, holes, summary
# For Rescue: char_rescue
```

### Complete Examples

#### Example 1: Basic CDD Processing
```bash
python domain_extract.py -c singles.cdd -o output/ -d results.csv
```

#### Example 2: Comprehensive Analysis with All Plots
```bash
python domain_extract.py -c singles.cdd -m 1 -p all -o analysis_results/ -d complete_data.csv
```

#### Example 3: Targeted Family Analysis
```bash
python domain_extract.py -c singles.cdd -t "2.A.1,3.A.10" -p "arch,summary" -o targeted_output/
```

#### Example 4: Rescue Data Processing
```bash
python domain_extract.py -r rescued/ -p char_rescue -o rescue_results/ -d rescue_data.csv
```

#### Example 5: Processing TCIDs from File
```bash
# Create a file with target TCIDs
echo -e "1.A.12\n2.A.1\n3.A.10" > target_families.txt

# Process the specified families
python domain_extract.py -c singles.cdd -t target_families.txt -p all -o targeted_analysis/
```

## Output Files

The tool generates several types of output:

### Generated Plots
Located in subdirectories under your output directory:
- `plots/general/` - General domain distribution plots
- `plots/char/` - Characteristic analysis plots  
- `plots/arch/` - Domain architecture visualizations
- `plots/holes/` - Hole detection and analysis plots
- `plots/summary/` - Summary statistics plots
- `plots/resc/` - Rescue-specific plots (when using rescue input)

### Data Files
- **CSV Export**: Structured data with columns for Accession, Length, Family, Subfamily, Domains, and Separators
- **Processing Logs**: Console output showing progress and execution time

## Input Data Requirements

### CDD File Format
- Standard CDD (Conserved Domain Database) format
- Contains domain hit information for protein sequences

### Rescue Directory Structure
- Directory containing rescue files
- Each file should contain domain information for protein families

### TCID Format
- Family format: `X.X.X` (e.g., `2.A.1`, `1.C.105`)
- Subfamily format: `X.X.X.X` (e.g., `2.A.1.1`, `3.A.10.11`)

## Troubleshooting

### Common Issues

1. **File Not Found Error**
   ```bash
   # Ensure your input file exists
   ls -la your_file.cdd
   ```

2. **Permission Issues**
   ```bash
   # Make sure output directory is writable
   chmod 755 output/
   ```

3. **Memory Issues with Large Files**
   - Process specific families instead of all families
   - Use the `-t` option to target specific TCIDs

4. **Missing Dependencies**
   ```bash
   # Reinstall requirements
   pip install -r requirements.txt --upgrade
   ```

### Performance Tips

- Use targeted TCID processing (`-t`) for faster execution on large datasets
- Specify output directory (`-o`) to organize results
- Use specific plot options instead of "all" for faster processing

## Examples of Advanced Usage

### Batch Processing Multiple Families
```bash
# Create a script for batch processing
for family in "1.A.12" "2.A.1" "3.A.10"; do
    python domain_extract.py -c singles.cdd -t "$family" -p arch -o "output_${family}/"
done
```

### Combining Results
```bash
# Process with data export and combine results
python domain_extract.py -c singles.cdd -t family_list.txt -d batch_results.csv -o batch_output/
```

## Project Structure

```
TCDB-Domain-Visualization/
├── domain_extract.py      # Main CLI tool
├── Family.py             # Family processing class
├── RescueFamily.py       # Rescue family processing class
├── requirements.txt      # Python dependencies
├── utility/              # Helper modules
│   ├── cdd_parser.py    # CDD file parsing
│   ├── rescue_parser.py # Rescue file parsing
│   ├── util.py          # General utilities
│   └── config.py        # Configuration settings
├── sequences/           # Input sequence files
├── output/              # Default output directory
├── plots/               # Generated visualization plots
└── rescued/             # Rescue data files
```

## Contributing

When contributing to this project:
1. Follow the existing code style
2. Add appropriate documentation for new features
3. Test with both CDD and rescue input formats
4. Update this README for any new CLI options

## License

[Add your license information here]

---

For additional help or questions, please refer to the documentation in the source code or create an issue in the project repository.
