# Building Your Own Pose Prediction Container
> This document details requirements and tips for writing a docking container for SAMPL-challenges. For Python template files that follow this guide, please see [SAMPL-containers/pose_prediction/templates/docking](https://github.com/samplchallenges/SAMPL-containers/tree/main/pose_prediction/templates/docker). For an example of a run-able docking main file, please see [SAMPL-containers/adv/main.py](https://github.com/samplchallenges/SAMPL-containers/blob/main/pose-prediction/examples/adv-docker/main.py). This guide is written under the assumption the reader has already gone through the [Docking Tutorial](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md). 

## Input Requirements
> Every container must be able to handle the following input arugments. These are the only inputs your container will be expected to handle. We recommend using [`Click`](https://pypi.org/project/click/) or [`argparse`](https://docs.python.org/3/library/argparse.html) to handle command line argument parsing.

### Required Inputs
* Your container will be run with all the following inputs:
  * `container-name:tag --receptor [file] --smiles [str] --hint [file] --hint-radius [float] --hint-molinfo [str] --output-dir [path]`
* `--receptor`: receptor `.pdb` file to dock the ligand into
  * Example: `--receptor data/receptor.pdb`
* `--smiles`: quoted SMILES string representing the ligand to dock (i.e. "CCC")
  * Example: `--smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O"`
* `--hint`: `.pdb` file with a hint ligand to denote the docking region
  * Example: `--hint data/hint.pdb`
* `--hint-radius`: numeric (float) Angstrom distance from the hint ligand (see above) to consider as the docking region
  * Example: `--hint-radius 6.0`
* `--hint-molinfo`: resname of the hint ligand used in the hint `.pdb` file
  * Example: `--hint-molinfo "E4Y"`
* `--output-dir`: directory to save final docking files (docked ligand and receptor files)
  * You will not need to handle determining the output directory input as the [`ever_given`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/ever_givenUsage.md) wrapper handles this for you. Please ensure that the required output files are saved to the directory path specified by the `output-dir` argument

### Optional Inputs
> The Optional Inputs will be most helpful to participants who are using proprietary licenses or files that cannot be uploaded to a public repository. It will also be helpful to participants who would like to explicilty separate their licenses from their container so their container cannot be run by just anyone. Please see the documentation on [License and Code Privacy](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/LicenseAndCodePrivacy.md) for more information.
* The [app.samplchallenges.org](https://app.samplchallenges.org/) Submission Form provides a "Special Arguments" section that allows you to specify command line arguments and a corresponding file to be passed to your container at runtime.  
* Any optional input arguments should be in the form of `--your-argument` with no capital letters. The argument should expect a file that you will upload as the input value.
* Your container will be run with the following inputs in the format:
  * `container-name:tag --receptor [file_path] --smiles [str] --hint [file_path] --hint-radius [float] --hint-molinfo [str] --output-dir [path] --your-argument [file_uploaded_by_you]`

## Output Requirements

**File Outputs**: Output the following files into the `output-dir`
* **docked ligand file**: a `.mol2`, `.pdb` or `.sdf` of the docked ligand
  * path_to_docked_ligand_file = `{output_dir}/{docked_ligand_file}`
* **receptor file**: a `.pdb` file of the receptor used or modified by your docking program, this is important for rmsd scoring purposes in case your complex changes frame of reference during the docking program
  * path_to_receptor_file = `{output_dir}/{receptor_file}`


**Printed Outputs**: Print the following to to `stdout`
* The LAST lines your container should output to `stdout` are below in the format `key value` where the keys are `docked_ligand` & `receptor` and the values are file paths. The key and value should be separated by a single space. You may print other outputs to `stdout` throughout your program, but these two lines must be the LAST lines printed by your program.
   ```
   docked_ligand {path_to_docked_ligand_file}
   receptor {path_to_receptor_file}
   ```

**Intentional No Prediction Output**: Output the following file and print the following to `stdout`
* If you are intentionally avoiding a prediction for a compound, please replace `{path_to_docked_ligand_file}` and `{path_to_receptor_file}` with a file saved to `output_dir` called `no_prediction` (see example below). If your output files are called `no_prediction`, any contents of the file will be ignored. 
   * path_to_no_prediction = `{output_dir}/no_prediction`
   ```
   docked_ligand {path_to_no_prediction}
   receptor {path_to_no_prediction}
   ```

## Program Logs
* Any output to `stdout` or `stderr` will be logged with timestamps associated with each output. These logs will be made accessible to you.
* Please print general logging info to `stdout` and error messages to `stderr` as is convention.
* As stated in [OutputRequirements](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/BuildYourOwnPosePredictionContainer.md#output-requirements), the last two lines of `stdout` output must be your two `key value` pairs. 


## A Note about Intermediate Files
* Please do not use `output_dir` to store your intermediate files
* You can store your intermediate files in a temporary directory or in any directory other than the `output_dir`, these files should die when your container finishes executing.


## Example Python Main Function Definition
> Every docking container you build for SAMPL challenges should include a main file with a main function that looks similar to the code block below. The following docking main template meets all input and output requirements mentioned above. 
```
import click
import os.path

@click.command()
@click.option("--receptor",required=True,type=click.Path(exists=True),help="path of receptor PDB to dock the ligand into")
@click.option("--smiles",required=True,help="SMILES str of ligand to be docked. quote to prevent CLI errors \"CCC\"")
@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint-molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint-radius",required=True,type=float,help="box size of the box to dock into")
@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")
def docking_main(receptor, smiles, hint, hint_molinfo, hint_radius, output_dir):
        ''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:        file    receptor PDB path to dock ligand into
                       smiles:          str     SMILES string of ligand to be docked 
                       hint:            file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo:    str     resname of the ligand used in the hint PDB
                       hint_radius:     float   radius around the hint ligand to consider in docking
                       output_dir:      str     output director for receptor and docked_ligand\
		       
            OUTPUTS:   prints           docked_ligand {path_to_docked_ligand_file}
	    				where the path_to_docked_ligand_file = output_dir/docked_ligand_file
                       prints           receptor {path_to_receptor_file}
		       			where the path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
                       writes file      docked ligand file as a .pdb .mol2 or .sdf this file must be saved on disk to the
		       			path stored in the output_dir argument
                       writes file      receptor prepped and used by program in docking as .pdb this file must be saved 
		       			on disk to the path stored in the output_dir argument
		       
	    RETURNS:   None
	    	       Please note that anything your function returns will be ignored by our automated scoring 
		       All outputs (docked_ligand and receptor files) MUST BE saved on disk to the specified
		       output_dir AND the absolute path must be printed out as specified in "OUTPUTS"
        '''
        
    
        # set the output file names / paths
        print("logging: setting filenames and paths")
        
        docked_ligand_file_name = ""
        receptor_file_name = ""
        path_to_docked_ligand_file = os.path.join(output_dir, docked_ligand_file_name)
        path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
       
        
        # YOUR DOCKING CODE GOES HERE
        
        # write out the docked ligand file to path_to_docked_ligand_file, your docked ligand file must be 
	# saved to the path specified by the output_dir parameter 
	
        # write out the receptor file to path_to_receptor_file, your receptor file must be 
        # saved to the path specified by the output_dir parameter         
        
        
        
        # Your final ligand and receptor files should be SAVED to the output_dir (specified as a parameter). 
        # Your main function should also PRINT out the 'key value' pairs 
        #	 * key: either the 'docked_ligand' or 'receptor'. 
        # 	 * value: the absolute file path to the file on disk. The absolute file paths to your files should be specified
        #	          as 'output_dir/your_file' where output_dir is the output_dir (specified as a parameter). Your final
        # 		  output files MUST BE saved to the output_dir, otherwise our automated scoring will not be able to 
        #		  find your files
        # 		       * path_to_docked_ligand_file = os.path.join(output_dir,"rec-dock.pdb") 
        # 		       * path_to_receptor_file = os.path.join(output_dir,"best_dock.pdb")
        
        
        print(f"docked_ligand {path_to_docked_ligand_file}")
        print(f"receptor {path_to_receptor_file}")
        
        
        # A NOTE ABOUT RETURN STATEMENTS
        # Anything your container returns will be ignored. Please make sure that any outputs follow the 'key value'
        # format mentioned in the previous comment
```

## Including your own Python Modules
If you modularize your code and include your own python modules, you will need to follow the steps below. For an example with using extra python modules beyond just main.py, please see [SAMPL-containers/docking/examples/adv-tutorial](https://github.com/samplchallenges/SAMPL-containers/tree/main/pose_prediction/examples/adv-docker).
1. Write your own python module(s)
2. Copy them into your Docking container using the `COPY` command in your Dockerfile or Singularity Definition File
    * Dockerfile:
    	```
	COPY main.py setup.py {your_python_module} ./
	```
    * Singularity Definition File:
    	```
	%files
	main.py {destination}
	setup.py {destination}
	{your_python_module} {destination}
3. Include your docking modules in the `py_modules` section of `setup.py`
    ```
    py_modules=[
       'main',
       '{your_python_module}',
    ]
    ```
    
## Including your main function as the ENTRYPOINT
If you use different naming conventions than those used in the template files for your main py file and main function, you will need to follow the steps below.
1. Write your own main py module and main function using your own naming conventions
2. Include your docking main in the `py_modules` section of `setup.py` 
    ```
    py_modules=[
       '{your_py_main}',
    ]
    ```
3. Alter the `entry_point` in `setup.py` to match you naming convention
    ```
    entry_points='''
        [console_scripts]
        {your_entrypoint_name}={your_py_main}:{your_main_function}
    '''
    ```
4. Copy the file into your Docking container
   * Docker: 
		```
		COPY {your_py_main} ./
		```
   * Singularity:
		```
		%files
		{your_py_main} {destination_file_path}
		```
5. Add your `entry_point` from step 3 to your Dockerfile or buildfile
   * Docker: 
		```
		ENTRYPOINT ['{your_entrypoint_name}']
		```
   * Singularity:
		```
		%runscript
		exec {your_entrypoint_name} $@
		```

## Tips for modifying the docking tutorial to fit your needs
> In some cases, the [miniconda3 docker image](https://hub.docker.com/r/continuumio/miniconda3) specified in the [tutorial `Dockerfile`](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#13-install-conda-environment-from-section-12-into-your-container) will not be compatible with the programs you your docking container will require. When this is the case, please try the following: 
* Go to [dockerhub](https://hub.docker.com/) and use the search bar to search for a container that meets your needs. 
	* For example, if I needed a container with a `gcc` compiler I would search for `gcc`, choose an image and locate the image name:
		![searchbar](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/images/dockerhub_search.png)
* Once you have the name of the image you will use as your base, let's call it `image-to-use`, change your Dockerfile or Singularity Definition File to build off of your `image-to-use`
	* Docker:
		```
		FROM image-to-use
		```
	* Singularity:
		```
		Bootstrap: docker
		From: image-to-use
		```
	
* If you still need `conda` for virtual environment management, we recommend installing `miniconda` by adding the installation steps into the Dockerfile or Singularity Definition File. An example is in the code blocks below:
	* Docker:
		```
		RUN wget \
		    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
		    && mkdir /root/.conda \
		    && bash Miniconda3-latest-Linux-x86_64.sh -b \
		    && rm -f Miniconda3-latest-Linux-x86_64.sh
		```
	* Singularity
		```
		%post
		wget \
		    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
		    && mkdir /root/.conda \
		    && bash Miniconda3-latest-Linux-x86_64.sh -b \
		    && rm -f Miniconda3-latest-Linux-x86_64.sh
		```
* Any additional steps to install your required programs should also be added to your `Dockerfile` for Singularity Definition File for more information, please see:
	* Docker: [Section 1.4](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#14-download-and-prepare-the-command-line-programs-autodock-vina-and-mgl-tools-executables-for-use-in-the-docking-container) and [Section 1.5](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#15-install-autodock-vina-and-mgl-tools-into-your-container) of the tutorial. 
	* Singularity: []() and []() of the tutorial


## Tips for integrating command line programs
* Some common command line programs (such as AutoDock Vina) already have docker containers made by other people or organizations. It may be worth it to search for pre-made docker containers to inherit from or build off of. (see [AutoDock Vina Docker](https://hub.docker.com/r/taccsciapps/autodock-vina))
* Some common command line programs may also have Python APIs (see [AutoDock Vina API](https://pypi.org/project/vina/)) 
* If the above bullets do not work, you can install the command line program into your container by copying the files into the container and running the installation steps in the Dockerfile or Singularity Definition File
    * Docker: [`SAMPL-league/docking/examples/adv-base/Dockerfile`](https://github.com/samplchallenges/SAMPL-containers/blob/main/pose_prediction/examples/adv-docker/Dockerfile)
    * Singularity: []()
* To run a command line program from within a Python module, consider using the [`subprocess`](https://docs.python.org/3/library/subprocess.html) library or [`os.system()`](https://docs.python.org/3/library/os.html?highlight=os%20system#os.system) or similar from the Python3 library
    * Please see [`SAMPL-league/docking/examples/adv/autodock.py`](https://github.com/samplchallenges/SAMPL-containers/blob/main/pose_prediction/examples/adv-docker/autodock.py)
