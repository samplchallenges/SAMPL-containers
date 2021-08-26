```
import click
from rdkit import Chem
from rdkit.Chem import Crippen

@click.command()
@click.option("--solute", default="CCCCCCCCO", help="solute SMILES string")
@click.option("--solventa", default="O", help="solventA SMILES string")
@click.option("--solventb", default="CCCCCCCCO", help="solventB SMILES string")
@click.option("--output-dir", help="Output Directory", type=click.Path(exists=True))
def get_logd(solute, solventa, solventb, output_dir):
    '''   INPUTS   solute      str   SMILES string representing the solute
                   solventa    str   SMILES string representing Solvent A
                   solventb    str   SMILES string representing Solvent 
                   output-dir  str   Directory to output files to (should not be used in this challenge)
    '''
    logP = None
    
    
    # YOUR LOGP CALCULATION CODE HERE
    

    print(f"LogD {logP}")
