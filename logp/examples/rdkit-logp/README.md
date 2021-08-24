# Description:
Calculates LogP of the solute SMILES string using `rdkit`'s LogP calculation

# Setup:
Not applicable


# Build:
1. `cd SAMPL-league/attic/predict-rdkitlogp`
2. `docker build -t rdlogp .`


# Run:
### Options
```
docker run -it --rm rdlogp --help
Usage: print-LogP [OPTIONS]

Options:
  --fuzz             Random change logP value by +/- 10%  [default: False]
  --smiles TEXT      ligand SMILES string
  --output-dir PATH  Output Directory
  --help             Show this message and exit.
```
### Run commands
`python ../../ever_given/run.py rdlogp --smiles <SMILES_str>`
* `python ../../ever_given/run.py rdlogp --smiles "CCCCc1ccccc1"`
