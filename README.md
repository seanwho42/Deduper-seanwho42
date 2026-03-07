# Deduper

Deduper is a PCR deduplication I was assigned to create as part of the Bi624/Bi625 course at UO KCGIP Bioinformatics and Genomics Masters Program. It supports error correction with a Hamming distance of 1 (further customization of this is a potential plan for the future). Only single-end reads are supported at this time.


## Installation

Deduper is installable through GitHub.

```bash
git clone https://github.com/seanwho42/Deduper-seanwho42.git
```

## Usage

Input sorted SAM file and text file containing known UMIs of a fixed length separated by new lines.

Outputs a copy of provided SAM file with PCR duplicates removed.

```
options:
  -h, --help        show this help message and exit
  -f, --file FILE   Sorted SAM file to have PCR duplicates removed from
  -o, --out OUT     Output path to sorted sam file without PCR duplicates
  -u, --umi UMI     File path for file containing the list of valid UMIs
  -c, --correction  Enable error correction to find the closeset fuzzy match (using Hamming distance <= 1)
```
