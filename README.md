# Mitochondrial-Genome-Rotation
A Python tool for rotating mitochondrial genomes based on the tRNA-Pro gene.


## Features

- Reads mitochondrial genome files in FASTA and GenBank formats.
- Detects the tRNA-Pro gene and determines the rotation position.
- Rotates the circular genome sequence.
- Updates gene and feature coordinates, including features that cross the sequence boundary.
- Validates sequence length and feature consistency.
- Exports the rotated genome in FASTA and GenBank formats.

## Technologies

- Python
- Biopython

## How to Run

1. Install Biopython:
   ```bash
   pip install biopython
   ```
2. Run the program:
   ```bash
   python main.py
   ```
3. Enter the FASTA and GenBank file paths when prompted.

## Output

The program generates two files:

- `rotated_genome.fasta` — Rotated mitochondrial genome sequence.
- `rotated_genome.gb` — Rotated genome with updated annotations.
