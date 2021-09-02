# Building Your Own LogD Container
> This document details requirements and tips for writing a LogD container for SAMPL-challenges.

## Input Requirements
Every container must be able to handle the following input flags. These are the only flags your container will be expected to handle. We typically use [`click`](https://click.palletsprojects.com/en/8.0.x/) to handle command line argument parsing, but feel free to use your preferred parser.
* `--solute`: smiles string representing the solute
* `--solventa`: smiles string representing solvent A 
* `--solventb`: smiles string representing solvent B

## Output Requirements
**Printed Outputs**: Print the following to to `stdout`
The LAST line your container should output is below in the format key value where the key is logd and the value is the float LogD value. The key and value should be separated by a single space. You may print other outputs throughout your program, but this line must be the LAST line printed by your program.
```
logd {logd_float_value}
```

## Program Logs
* Any output to `stdout` or `stderr` will be logged with timestamps associated with each output. These logs will be made accessible to you.
* Please print general logging info to `stdout` and error messages to `stderr` as is convention.
* Feel free to print to `stdout` as needed, but as stated in [OutputRequirements](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/DockingContainerRequirements.md#output-requirements), the last line of output must be your `key value` pair. 

## Example Python Main Function Definition
```
import click

@click.command()
@click.option("--solute", help="smiles string representing the solute")
@click.option("--solventa", help="smiles string representing solvent A")
@click.option("--solventb", help="smiles string representing solvent B")

def logd_main(solute, solventa, solventb):
        ''' calculates LogD the given smiles strings for solute, solvent a and solvent b
            INPUTS:    solute           str    smiles string representing the solute
                       solventa         str    smiles string representing solvent A 
                       solventb         str    smiles string representing solvent B
            OUTPUTS:   prints           logd {logd_float_value}
        '''
        LogD = None
        print("logging: calculating LogD")
        
        # YOUR LOGD CODE GOES HERE
        
        
        
        print("logging: finished calculating LogD")
        
        # print out the key value LogD pair. Please note that this is the last thing the program outputs
        print(f"logd {LogD}")
```
