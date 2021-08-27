# LogD Container Requirements
> This document details the input and output requirements for building and running a LogD container for SAMPL-challenges. 

## Input Requirements
* `--solute`: SMILES string representing the solute
* `--solventa`: SMILES string representing solvent A which forms part of the partition
* `--solventb`: SMILES string representing solvent B which forms part of the partition
* `--output-dir`: output directory path for any output files, this parameter will not be used for the LogD challenge

## Output Requirements
* **Printed Output**: print the following to `stdout`
    ```
    LogD {LogD_float}
    ```
    * Your container should output above in the format key value where the key is `LogD` and the value is the caluclated LogD value as a float. The key and value should be separated by a single space
    * These are the only two outputs that should be printed to `stdout`. Please print any extraneous error messages to stderr so output parsing is not compromised
    * If you are purposely avoiding outputtting a prediction for a compound, please replace `{LogD_float}` with `no_prediction` (see example below)
        ```
        LogD no_prediction
        ```


## Example LogD Main Function
> Every LogD container you build for SAMPL challenges should include a main file that looks something code block below. The following docking main template meets all input and output requirements mentioned above.
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
    logD = None
    
    
    # YOUR LOGP CALCULATION CODE HERE
    

    print(f"LogD {logD}")
    
```


## Including your own Python Modules
If you modularize your code and include your own python modules, you will need to do the following:
1. Write your own python module(s)
2. Copy them into your LogP container using the `COPY` command in your Dockerfile
    * `COPY main.py setup.py <your_python_module>.py ./`
3. Include your docking modules in the `py_modules` section of `setup.py`
    ```
    py_modules=[
       'main',
       '{your_python_module}',
    ]
    ```
