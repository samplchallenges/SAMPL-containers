import click

@click.command()
@click.option("--solute", help="smiles string representing the solute")
@click.option("--solventa", help="smiles string representing solvent A")
@click.option("--solventb", help="smiles string representing solvent B")
@click.option("--output-dir", help="disregard for this container")

def logd_main(solute, solventa, solventb, output_dir):
        ''' calculates LogD the given smiles strings for solute, solvent a and solvent b
            INPUTS:    solute       str    smiles string representing the solute
                       solventa     str    smiles string representing solvent A 
                       solventb     str    smiles string representing solvent B
                       output_dir   str    disregard for this container
            OUTPUTS:   prints       logd {logd_float_value}
        '''
        LogD = None
        print("logging: calculating LogD")
        
        # YOUR LOGD CODE GOES HERE
        
        
        
        print("logging: finished calculating LogD")
        
        # print out the key value LogD pair. Please note that this is the last thing the program outputs
        print(f"logd {LogD}")
