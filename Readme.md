# ClaPNAC: A Classifier of Protein–Nucleic Acid Contacts

ClaPNAC is a geometric classifier for identifying and classifying contacts between amino acid (AA) and nucleic acid (NA) residues in three-dimensional molecular structures.

The tool annotates protein–nucleic acid complexes with pairwise contacts and assigns interaction classes based on the spatial arrangement of amino acid and nucleotide atoms. Unlike methods that focus primarily on specific chemical interactions, ClaPNAC classifies contacts using the geometry of interacting residues.

ClaPNAC can be applied to RNA–protein and DNA–protein complexes, as well as isolated nucleotides, nucleosides, and amino acids.

**Version:** 1.00

**Preprint:** Nikolaev G, Baulin EF, Bujnicki JM. *ClaPNAC: a Classifier of Protein–Nucleic Acid Contacts*. bioRxiv (2025). https://doi.org/10.1101/2025.05.30.657037

**Repository:** https://github.com/GrigorijNik/ClaPNAC

## Table of Contents

- [Features](#features)
- [Interaction Classification](#interaction-classification)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#input-data)
- [Usage](#usage)
- [Command-Line Arguments](#command-line-arguments)
- [Examples](#examples)
- [Output](#output)
- [Comparing Results](#comparing-results)
- [Scientific Background](#scientific-background)
- [Citation](#citation)
- [License](#license)

## Features

ClaPNAC provides the following functionality:

- Identification and classification of amino acid–nucleotide contacts in three-dimensional molecular structures.
- Classification of contacts into seven interaction classes.
- Processing of PDB and mmCIF structure files.
- Analysis using full-atom (FA) or coarse-grained (CGR) representations.
- Selection of a specific structural model for analysis.
- Filtering of predicted contacts using a configurable score threshold.
- Control of the spatial range used for contact extraction.
- Optional sequence-based output highlighting residues involved in interactions.
- Optional reporting of only the highest-scoring interaction for each residue pair.
- Comparison of previously generated results using set difference or union operations.

## Interaction Classification

ClaPNAC distinguishes seven classes of amino acid–nucleotide interactions:

| Class | Description |
|---|---|
| `P` | Interaction with the nucleotide phosphate group |
| `R` | Interaction with the nucleotide ribose |
| `W` | Interaction with the Watson–Crick edge of the nucleobase |
| `H` | Interaction with the Hoogsteen edge of the nucleobase |
| `S` | Interaction with the sugar edge of the nucleobase |
| `+` | Stacking interaction involving one face of the nucleobase |
| `-` | Stacking interaction involving the opposite face of the nucleobase |

Each classified interaction is associated with a score between 0 and 1. The score reflects the geometric compatibility of the input residue pair with the corresponding interaction class.

ClaPNAC uses a database of structurally characterized amino acid–nucleotide doublets to evaluate potential interactions.

The output can contain multiple interaction classes for a single amino acid–nucleotide residue pair.

### Representative Structures for Interaction Classes

Representative examples of amino acid–nucleotide interactions are provided in the `classes_PDB/` directory.

The directory is organized by amino acid residue type. Each amino acid folder contains subdirectories with representative structures classified into different interaction classes.

For example:

```text
classes_PDB/
├── ALA/
│   ├── select_class_A-
│   ├── select_class_A+
│   ├── select_class_AH
│   ├── select_class_AS
│   ├── select_class_AW
│   ├── select_class_C-
│   ├── select_class_C+
│   ├── select_class_CH
│   ├── select_class_CS
│   ├── select_class_CW
│   ├── select_class_G-
│   ├── select_class_G+
│   ├── select_class_GH
│   ├── select_class_GS
│   ├── select_class_GW
│   ├── select_class_P
│   ├── select_class_R
│   ├── select_class_U-
│   ├── select_class_U+
│   ├── select_class_UH
│   ├── select_class_US
│   └── select_class_UW
├── ARG/
├── ASN/
├── ASP/
└── ...
```

The same organization is used for other amino acid residue types, including ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, and VAL.

The subdirectory names indicate the interaction class and, where applicable, the nucleobase identity. For example, `select_class_AW` represents an interaction involving adenine and the Watson–Crick edge, whereas `select_class_P` and `select_class_R` correspond to phosphate and ribose interactions, respectively.

These representative structures can be used to explore the geometric features underlying each interaction class and to understand how ClaPNAC classifies amino acid–nucleotide contacts.

## Requirements

ClaPNAC is implemented in Python.

The dependencies include:

- Python 3.8
- NumPy
- pandas
- SciPy
- Biopython
- argparse
- seaborn
- matplotlib
- alive-progress

Some dependencies, including `argparse`, are part of the Python standard library.

The exact dependency versions should be selected according to the versions supported by the release.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/GrigorijNik/ClaPNAC.git
cd ClaPNAC
git checkout v1.00
```

### 2. Create a Python environment

Using Conda:

```bash
conda create -n clapnac python=3.8
conda activate clapnac
```

### 3. Install dependencies

Install the required Python packages:

```bash
pip install numpy pandas scipy biopython seaborn matplotlib alive-progress
```

## Input Data

ClaPNAC accepts molecular structures in PDB (`.pdb`) and mmCIF (`.cif`) formats.

Create a directory containing the structures to be analyzed. The program recursively searches the specified directory and its subdirectories for files with these extensions.

For example:

```text
input/
├── complex_1.pdb
├── complex_2.cif
└── additional_structures/
    ├── complex_3.pdb
    └── complex_4.cif
```

The input structures should contain the amino acid and nucleic acid residues of interest.

When a structure contains multiple models, the desired model can be specified using the `-m` option.

## Usage

ClaPNAC is executed from the command line.

The general command format is:

```bash
python ClaPNAC.py \
    -i <input_directory> \
    -t <FA_or_CGR> \
    -o <score_threshold> \
    -r <contact_range>
```

The remaining arguments specify the input directory, atom representation, minimum output score, and contact range.

Optional arguments enable sequence output, model selection, unique-interaction reporting, and comparison of previously generated results.

## Command-Line Arguments

| Argument | Description | Default |
|---|---|---|
| `-i`, `--in_dir` | Directory containing input PDB and/or mmCIF files | Required for structure analysis |
| `-t`, `--type_in` | Atom representation: `FA` (full atom) or `CGR` (coarse-grained) | `FA` |
| `-o`, `--output` | Minimum score threshold for reporting interactions; a floating-point value between 0 and 1 | `1.0` |
| `-r`, `--range_in` | Spatial range used by ContExt for contact extraction, in Å | `3.9` |
| `-m`, `--model` | Structural model identifier, for example `#1` | `#1` |
| `-s`, `--seq` | Enables sequence-based output, marking interacting residues with `*` and other residues with `.` | Disabled |
| `-u`, `--unique` | Reports only the highest-scoring interaction for each residue pair | Disabled |
| `-c`, `--comparison` | Compares previously generated results using `D` (difference) or `U` (union) | Disabled |

### Atom representation (`-t`)

Two representations are supported:

- `FA`: full-atom representation.
- `CGR`: coarse-grained representation.

The selected representation determines how the structure is represented during contact classification.
Atoms used for the coarse-grained representation could be found in classes/CGR.py.

### Score threshold (`-o`)

The score threshold determines which classified interactions are included in the output.

For example:

- `-o 1.0`: reports interactions with a score of 1.0.
- `-o 0.8`: reports interactions with a score of at least 0.8.
- `-o 0.5`: reports interactions with a score of at least 0.5.

Lowering the threshold allows more geometrically compatible residue pairs to be reported.

### Contact range (`-r`)

The contact range specifies the spatial range used by ContExt during contact extraction. The default value is `3.9` Å.

### Structural model (`-m`)

The model option specifies which structural model should be analyzed. The default identifier is `#1`.

For structures containing multiple models, specify the identifier corresponding to the desired model.

### Sequence output (`-s`)

This option enables sequence-based output.

Residues involved in interactions are marked with `*`, while residues without reported interactions are marked with `.`.

For example:

GGACUAGCGGAGGCUAGUCC
..***..*******......

QKGNFRNQRKTVKCFNCGKEGHIAKNCRAPRKKGCWKCGKEGHQMKDCTERQAN
...**.*...*..*.*.........*.....*..*.*.................

### Unique interactions (`-u`)

By default, a residue pair may be associated with multiple interaction classes.

When `-u` is enabled, only the highest-scoring interaction (above the score threshold) is retained for each residue pair.

This option is useful when a simplified representation of the interaction network is required.

## Examples

### 1. Analyze structures using the default settings

```bash
python ClaPNAC.py -i ./input
```

The program processes PDB and mmCIF files found in the input directory and its subdirectories using the default full-atom representation, a score threshold of 1.0, and a contact range of 3.9 Å.

### 2. Analyze structures using the full-atom representation

```bash
python ClaPNAC.py \
    -i ./input \
    -t FA \
    -o 1.0 \
    -r 3.9
```

### 3. Analyze structures using the coarse-grained representation

```bash
python ClaPNAC.py \
    -i ./input \
    -t CGR \
    -o 0.8 \
    -r 3.9
```

### 4. Include interactions with scores of at least 0.5

```bash
python ClaPNAC.py \
    -i ./input \
    -t FA \
    -o 0.5
```

### 5. Generate sequence-based output

```bash
python ClaPNAC.py \
    -i ./input \
    -s
```

### 6. Report only the highest-scoring interaction for each residue pair

```bash
python ClaPNAC.py \
    -i ./input \
    -u
```

### 7. Combine sequence output and unique-interaction reporting

```bash
python ClaPNAC.py \
    -i ./input \
    -t FA \
    -o 0.5 \
    -s \
    -u
```

### 8. Select a structural model

```bash
python ClaPNAC.py \
    -i ./input \
    -m "#1"
```

The example selects model `#1`. Use the model identifier supported by the input structure and the implementation.

## Output

For each processed structure, ClaPNAC generates CSV output containing the classified amino acid–nucleotide contacts.

The output includes information describing the interacting residues, classification scores, and interaction types.

The principal output fields are:

| Field | Description |
|---|---|
| `NA_id` | Identifier of the nucleic acid residue |
| `P_id` | Identifier of the protein residue |
| `score` | Score assigned to the interaction class |
| `type` | Interaction class (`P`, `R`, `W`, `H`, `S`, `+`, or `-`) |
| `N-AA` | Nucleotide–amino acid residue pair |

An illustrative record is:

```text
NA_id   P_id   score   type   N-AA
B-203   A-4    0.63    P      A-GLY
```

This record indicates that the `A 203` from chain `B` and the `GLY 4` from chain `A` have a phosphate-class interaction with a score of 0.63.

The output filename encodes information about the input structure and analysis settings. Consult the generated filenames and the corresponding output files when processing multiple structures.

**Important:** The score is a geometric classification score. It should not be interpreted directly as a binding affinity or a calibrated probability of biological interaction.

## Comparing Results

ClaPNAC includes a comparison mode for comparing previously generated results.

The `-c` option accepts two values:

- `D`: difference.
- `U`: union.

The directory specified by `-i` is passed to the comparison utility.

### Difference

To calculate the difference between results:

```bash
python ClaPNAC.py \
    -i ./comparison_results \
    -c D
```

### Union

To combine results using the union operation:

```bash
python ClaPNAC.py \
    -i ./comparison_results \
    -c U
```

The comparison directory should contain the CSV files expected by the comparison utility.

## Scientific Background

Protein–nucleic acid interactions are essential for transcription, translation, DNA replication, genome maintenance, and the regulation of gene expression.

Understanding these interactions requires detailed characterization of contacts between individual amino acid and nucleotide residues.

ClaPNAC extends the geometric classification approach developed for ClaRNA to amino acid–nucleotide interactions. It uses structurally characterized residue pairs to classify contacts according to their spatial arrangement.

The method can be used to:

- Characterize protein–RNA and protein–DNA interfaces.
- Compare interaction patterns across related molecular structures.
- Identify changes in residue-level contacts between structurally similar complexes.
- Analyze interaction preferences in protein and nucleic acid families.
- Support the evaluation and interpretation of experimentally determined or computationally predicted molecular structures.

For the methodology, benchmark datasets, and scientific evaluation, please consult the preprint.

## Citation

If you use ClaPNAC in your research, please cite:

Nikolaev G, Baulin EF, Bujnicki JM. *ClaPNAC: a Classifier of Protein–Nucleic Acid Contacts*. bioRxiv, 2025.

DOI: https://doi.org/10.1101/2025.05.30.657037

Please also cite the software release when appropriate:

ClaPNAC, version 1.00. GitHub: https://github.com/GrigorijNik/ClaPNAC/tree/v1.00

## License

ClaPNAC is distributed under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

You are free to share and adapt the material for any purpose, including commercial use, provided that appropriate credit is given, a link to the license is provided, and any changes are indicated.

For the full license terms, see the Creative Commons Attribution 4.0 International license.

## Contact

For questions, bug reports, or suggestions, please use the GitHub repository:

https://github.com/GrigorijNik/ClaPNAC