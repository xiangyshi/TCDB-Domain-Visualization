# Input Data Requirements

The model expects input data with the following columns:

- **Protein**: Protein accession number
  - Example: `P0A334`

- **Protein_length**: Total length of the protein sequence
  - Example: `100`


- **Family**: Three-digit
  - Format: X.X.X (where X represents numbers)
  - Example: `2.A.1`, `2.A.11`


- **Subfamily**: Four-digit
  - Format: X.X.X.X (where X represents numbers)
  - Example: `2.A.1.1`, `3.A.10.11`



- **Domains**: List of tuples containing domain details
  - Format: List[(domain_accession, start_position, end_position, bit_score)]
  - Each tuple contains:
    - Domain accession
    - Start position
    - End position
    - Bit score
  - Example: `[(A, 1, 10, 8.71e-22), (B, 11, 20, 8.71e-22)]`

- **Motifs**: List of tuples containing motif details
  - Format: List[(motif_accession, start_position, end_position)]
  - Each tuple contains:
    - Motif accession
    - Start position
    - End position
  - Example: `[(A, 1, 10), (B, 11, 20)]`



- **Separator**: List of tuples describing spaces between domains in that protein
  - Format: List[(domain_transition, start_position, end_position)]
  - Each tuple contains:
    - Space information: String that describes from which domain to which domain (format: "X to Y")
    - Start position
    - End position
  - Example: `[("A to C", 78, 108), ("C to F", 179, 275)]`





















