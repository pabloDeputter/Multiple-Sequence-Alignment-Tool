## Multiple Sequence Alignment Tool

### Overview

This project implements a generalized dynamic programming algorithm for multiple sequence alignment (MSA), based on the Needleman-Wunsch and Smith-Waterman algorithms for pairwise alignment. The tool supports both global and local alignment of multiple DNA or protein sequences.

### Usage Instructions

#### Changing Parameters

Parameters such as scoring weights (match, mismatch, indel penalties) can be adjusted in the `data/settings.json` file.

#### Running the CLI

The CLI can be run with the default settings by executing:
```bash
python main.py -i data/input/input.fasta -o data/output/output.txt -a global
```

To use custom settings, ensure the `data/settings.json` file is updated accordingly before running the CLI:
```bash
python main.py -i data/input/input.fasta -o data/output/output.txt -s data/settings.json -a global
```

#### Input and Output Files

The input FASTA files should be placed in the `data/input` directory. The aligned output will be saved in the `data/output` directory. An example input file might look like:
```
>unknown_J_region_1
GYSSASKIIFGSGTRLSIRP
>unknown_J_region_2
NTEAFFGQGTRLTVV
>unknown_J_region_3
NYGYTFGSGTRLTVV
```
The output will be formatted as:
```
unknown_J_region_1: GYSSASKIIFGSGTRLSIRP
unknown_J_region_2: NTE.AF...FGQGTRLTVV.
unknown_J_region_3: NYG.YT...FGSGTRLTVV.
Alignment score: 45
```

#### Running Tests

The tests for the MSA algorithm can be run using:
```bash
python -m unittest tests/testMSA.py
```
This will execute all test cases and verify the correctness of the implementation against the expected outputs.

### Project Structure

- **data**: Contains input data files, output results, and the configuration file (`data/settings.json`) for scoring parameters (match, mismatch, gap penalties).
- **report.pdf**: The main report document.
- **requirements.txt**: Lists project dependencies.
- **src**: Holds the source code for the project, including the main script (`src/msa.py`) and other utility functions.
- **tests**: Contains the test cases for the project.
  - **data/expected**: Contains the expected output for the tests.
  - **data/input**: Contains the input data for the tests.
  - **test_msa.py**: Includes the tests for the MSA algorithm.
- **main.py**: The primary script for running the CLI. This file can be modified to run the MSA algorithm without the CLI by (un)commenting appropriate sections of code.

### Methodology

The generalized dynamic programming algorithm extends the Needleman-Wunsch and Smith-Waterman algorithms for pairwise alignment to k sequences, using a k-dimensional scoring matrix. The matrix is initialized with gap penalties for global alignment or zeros for local alignment. The scoring scheme considers all possible pairwise combinations of sequences. The matrix is filled iteratively, and the optimal alignment is reconstructed through a traceback procedure.

### Algorithmic Complexity

- **Space Complexity**: \( O(n^k) \)
- **Time Complexity**: \( O((2n)^k \cdot k^2) \)

### Application

The tool has been applied to biological problems, such as aligning T cell receptor sequences to identify conserved regions, demonstrating its practical utility in bioinformatics.