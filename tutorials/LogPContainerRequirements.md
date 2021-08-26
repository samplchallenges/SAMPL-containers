# LogD Container Requirements
> This document details the input and output requirements for building and running a LogD container for SAMPL-challenges.

## Input Requirements
* `--solute`: SMILES string representing the solute
* `--solventa`: SMILES string representing solvent A which forms part of the partition
* `--solventb`: SMILES string representing solvent B which forms part of the partition
* `--output-dir`: output directory path for any output files, this parameter will not be used for the LogD challenge

## Output Requirements


## Example LogD Main Function
```
import click

@click.command()
@click.option("--solute", default="CCCCCCCCO", help="solute SMILES string")
@click.option("--solventa", default="O", help="solventA SMILES string")
@click.option("--solventb", default="CCCCCCCCO", help="solventB SMILES string")
@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))
def logd_main(solute, solventa, solventb, output_dir):
    '''   INPUTS   solute      str   SMILES string representing the solute
                   solventa    str   SMILES string representing Solvent A
                   solventb    str   SMILES string representing Solvent 
                   output-dir  str   Directory to output files to (should not be used in this challenge)
    '''
    logP = None
    
    
    # YOUR LOGP CALCULATION CODE HERE
    

    print(f"LogD {logP}")
    
```
