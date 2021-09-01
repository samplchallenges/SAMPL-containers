# Building Your Own Docking Containter
> This document details the requirements for writing a docking container for SAMPL-challenges. For python template files that follow this guide, please see [SAMPL-containers/tutorials/template](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials/template). For a full example of a docking main file, please see [SAMPL-containers/adv/main.py](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/adv/main.py).

## Input Requirements
> Every container must be able to handle the following input flags. These are the only flags your container will be expected to handle. We typically use [`click`](https://click.palletsprojects.com/en/8.0.x/) to handle command line argument parsing, but feel free to use or build your own parser. 
* `--receptor`: receptor `.pdb` file to dock the ligand into
  * Example: `--receptor data/receptor.pdb`
* `--smiles`: quoted SMILES string representing the ligand to dock (i.e. "CCC")
  * Example: `--smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O"`
* `--hint`: `.pdb` file with a hint ligand to denote the docking region
  * Example: `--hint data/hint.pdb`
* `--hint_radius`: numeric (float) Angstrom distance from the hint ligand (see above) to consider as the docking region
  * Example: `--hint_radius 6.0`
* `--hint_molinfo`: resname of the hint ligand used in the hint `.pdb` file
  * Example: `--hint_molinfo "E51"`
* `--output-dir`: directory to save final docking files (docked ligand and receptor files)
  * You will not need to handle determining the output directory input as the [`ever_given`](https://github.com/samplchallenges/SAMPL-containers/blob/megosato-patch-1-1/tutorials/ever_givenUsage.md#ever_given-usage-notes) wrapper handles this for you. Please ensure that the required output files are saved to the output-dir directory

## Output Requirements

**File Outputs**: Output the following files into the `output-dir`
* **docked ligand file**: a `.mol2`, `.pdb` or `.sdf` of the docked ligand
  * path_to_docked_ligand_file = `{output-dir}/{docked_ligand_file}`
* **receptor file**: a `.pdb` file of the receptor used or modified by your docking program, this is important for rmsd scoring purposes in case your complex changes frame of reference
  * path_to_receptor_file = `{output-dir}/{receptor_file}`


**Printed Outputs**: Print the following to to `stdout`
* The LAST lines your container should output are below in the format `key value` where the keys are `docked_ligand`/`receptor` and the values are file paths. The key and value should be separated by a single space. You may print other outputs throughout your program, but these two lines must be the LAST lines printed by your program.
 ```
 docked_ligand {path_to_docked_ligand_file}
 receptor {path_to_receptor_file}
 ```
* These are the only two outputs that should be printed to `stdout`. Please print any extraneous error messages to `stderr` so output parsing is not compromised
* If you are purposely avoiding outputting a prediction for a compound, please replace `{path_to_docked_ligand_file}` and `{path_to_receptor_file}` with `no_prediction` (see example below)
   ```
   docked_ligand no_prediction
   receptor no_prediction
   ```

## Program Logs
* Any output to `stdout` or `stderr` will be logged with timestamps associated with each output. These logs will be made accessible to you.
* Feel free to print to `stdout` as needed, but as stated in [OutputRequirements](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/DockingContainerRequirements.md#output-requirements), the last two lines of output must be your two `key value` pairs. 


## Example Python Main Function Definition
> Every docking container you build for SAMPL challenges should include a main file that looks something code block below. The following docking main template meets all input and output requirements mentioned above. 
```
import click
import os.path

@click.command()
@click.option("--receptor", required=True, type=click.Path(exists=True), help="path of receptor PDB to dock the ligand into")
@click.option("--smiles", required=False, help="file with SMILES strings of ligands to be docked")
@click.option("--smiles_arg", required=False, help="SMILES string of a ligand within quotes to avoid command line parsing errors (i.e. \"CCC\")")

@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint_molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint_radius",required=True,type=float,help="box size of the box to dock into")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

def docking_main(receptor, smiles, smiles_argument, hint, hint_molinfo, hint_radius, output_dir):
        ''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:        file    receptor PDB path to dock ligand into
                       smiles:          str     SMILES string of ligand to be docked 
                       hint:            file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo:    str     resname of the ligand used in the hint PDB
                       hint_radius:     float   radius around the hint ligand to consider in docking
                       output_dir:      str     output director for receptor and docked_ligand
            OUTPUTS:   prints           docked_ligand {path_to_docked_ligand_file}
                       prints           receptor {path_to_receptor_file}
                       writes file(s)   docked ligand file as a .pdb .mol2 or .sdf
                       writes file(s)   receptor prepped and used by program in docking as .pdb
        '''
        
    
        # set the output file names / paths
        print("setting filenames and paths")
        docked_ligand_file_name = ""
        receptor_file_name = ""
        path_to_docked_ligand_file = os.path.join(output_dir, docked_ligand_file_name)
        path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
        
        
        
        # YOUR DOCKING CODE GOES HERE
        
        
        # write out the docked ligand file to path_to_docked_ligand_file
        print("writing docked ligand")
        
        # write out the receptor file to path_to_receptor_file
        print("writing prepared receptor file")
        
        
        # print out the key value pairs 
        #    * where the keys are docked_ligand and receptor
        #    * where the values are the file paths of the docked_ligand and receptor
        
        print(f"docked_ligand {path_to_docked_ligand_file}")
        print(f"receptor {path_to_receptor_file}")
```

## Including your own Python Modules
If you modularize your code and include your own python modules, you will need to follow the steps below. For an example with using extra python modules beyond just main.py, please see [SAMPL-containers/adv](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/adv).
1. Write your own python module(s)
2. Copy them into your Docking container using the `COPY` command in your Dockerfile
    * `COPY main.py setup.py <your_python_module>.py ./`
3. Include your docking modules in the `py_modules` section of `setup.py`
    ```
    py_modules=[
       'main',
       '{your_python_module}',
    ]
    ```
    
## Tips for integrating command line programs
* Some common command line programs (such as AutoDock Vina) already have docker containers made by other people or organizations. It may be worth it to search for pre-made docker containers to inherit from or build off of. (see [AutoDock Vina Docker](https://hub.docker.com/r/taccsciapps/autodock-vina))
* Some common command line programs may also have Python API's (see [AutoDock Vina API](https://pypi.org/project/vina/)) 
* If the above bullet doesn't work, you can install the command line program into your container by copying the files into the container and running the installation steps in the Dockerfile
    * Please see [`SAMPL-league/examples/adv-base/Dockerfile`](https://github.com/samplchallenges/SAMPL-league/blob/84ec83f00a637f9c79b6d2e3a1a336ea91837b7a/examples/adv-base/Dockerfile#L14)
* To run a command line program from within a Python module, consider using [os.system()](https://docs.python.org/3/library/os.html?highlight=os%20system#os.system) or similar from the Python3 library
    * Please see [`SAMPL-league/examples/adv/autodock.py`](https://github.com/samplchallenges/SAMPL-league/blob/84ec83f00a637f9c79b6d2e3a1a336ea91837b7a/examples/adv/autodock.py#L166)
