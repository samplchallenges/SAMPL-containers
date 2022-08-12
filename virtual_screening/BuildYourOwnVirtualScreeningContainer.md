# Building Your Own Virtual Screening Container
> This document details requirements and tips for writing a docking container for SAMPL-challenges. For Python template files that follow this guide, please see [SAMPL-containers/virtual_screening/templates/docker](https://github.com/samplchallenges/SAMPL-containers/tree/main/virtual_screening/templates/docker). For an example of a run-able virtual screening main file, please see [SAMPL-containers/virtual_screening/examples/adv-screen-docker/main.py](https://github.com/samplchallenges/SAMPL-containers/blob/main/virtual_screening/examples/adv-tutorial/main.py). This guide is written under the assumption the reader has already gone through the [Docking Tutorial](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md). 
> 
> 
> **IMPORTANT: Your container should be failsafe. If your program(s) fail on a specific output, you should still present a result/prediction. Perhaps you can predict a default value of 0 (non-binder) for whether the compound binds, and a value of 0 for the score. For ideas on how to make your container failsafe, please see [WritingFailsafeContainers.md]().**

## Input Requirements
> Every container must be able to handle the following input arugments. These are the only inputs your container will be expected to handle. We recommend using [`Click`](https://pypi.org/project/click/) or [`argparse`](https://docs.python.org/3/library/argparse.html) to handle command line argument parsing.

### Required Inputs
Your container will be run with all the following inputs:
* `container_name:tag [--batch] --receptor [file] --smiles [str/file] --c_x [float] --c_y [float] --c_z [float] --sz_x [float] --sz_y [float] --sz_z [float] `
* `--batch`: a flag that when present tells the container to expect batched inputs, when missing tells the container to expect single inputs. 
	* **Batched inputs**: in a single container run, the container will be expected to screen multiple ligands (passed as a text file of smiles) for whether they bind to the given protein. Batched inputs will be used to improve run time. Using batched inputs should only require the container to be loaded once per batch and the protein to be prepared once per batch.
	* **Un-batched inputs**: in a single container run, the container will be expected to screen a single ligands (passed as smiles string via command line) for whether it binds to the given protein.
* `--receptor`: receptor `.pdb` file to screen the ligand against
	* Example: `--receptor data/receptor.pdb`
* `--smiles`: 
	* **Batched inputs**: csv file in the format "id,name,value" where id=molid, name=molname, and value=smiles string
		```
		id,name,value
		0,SAMPL9-000,c1ccccc1
		1,SAMPL9-001,CCCCC
		```
	* **Un-batched inputs**: quoted SMILES string representing the ligand to dock (i.e. "CCC")
* `--c_x`: float value, center of docking box x coordinate
* `--c_y`: float value, center of docking box y coordinate
* `--c_z`: float value, center of docking box z coordinate
* `--sz_x`: float value, size of docking box in x direction
* `--sz_y`: float value, size of docking box in y direction
* `--sz_z`: float value, size of docking box in z direction

### Optional Inputs
> The Optional Inputs will be most helpful to participants who are using proprietary licenses or files that cannot be uploaded to a public repository. It will also be helpful to participants who would like to explicilty separate their licenses from their container so their container cannot be run by just anyone. Please see the documentation on [License and Code Privacy](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/LicenseAndCodePrivacy.md) for more information.
* The [app.samplchallenges.org](https://app.samplchallenges.org/) Submission Form provides a "Special Arguments" section that allows you to specify command line arguments and a corresponding file to be passed to your container at runtime.  
* Any optional input arguments should be in the form of `--your-argument` with no capital letters. The argument should expect a file that you will upload as the input value.
* Your container will be run with the following inputs in the format:
	* `container-name:tag --receptor [file_path] --smiles [str] --c_x [float] --c_y [float] --c_z [float] --sz_x [float] --sz_y [float] --sz_z [float] --output-dir [path] --your-argument [file_uploaded_by_you]`


## Output Requirements

**File Outputs**: 
* **Batched**: Output the following files by saving them into the `output-dir`
	* **scores csv file**: csv file in the format "id,name,value" where id=molid, name=molname, and value=docking score. 
		```
		id,name,value
		0,SAMPL9-000,-8.789
		1,SAMPL9-001,-5.349
		```
	* **binds csv file**: csv file in the format "id,name,value" where id=molid, name=molname, and value=whether the ligand binds (0 for non-binder, 1 for binder)
		```
		id,name,value
		0,SAMPL9-000,1
		1,SAMPL9-001,0
		```
* **Un-batched**: No file outputs

**Printed Outputs**:
* The LAST lines your container should output to `stdout` are below in the format `key value` where the keys are `docking_score` & `compound_binds` and the values are file paths. The key and value should be separated by a single space. You may print other outputs to `stdout` throughout your program, but these two lines must be the LAST lines printed by your program.
* **Batched**: print out each key with the file path to the csv 
	```
	docking_score {path_to_score_file_csv}
	compound_binds {path_to_binds_file_csv}
	```
* **Un-batched**: print out each key with the value
	```
	docking_score {float_docking_score}
	compound_binds {1_if_ligand_binds_0_if_nonbinder}
	```

**Intentional No Prediction Output**
* Intentional no predictions will not be permitted in this challenge. If you are unable to make a prediciton, please still output a score and a default prediction of whether the compound binds or not (1 or 0 respectively). Based on the makeup of this library, we recommend predicting the compound does not bind (value = 0) as a default.

## Running Your Container
If you haven't already, create a new virtual environment and install our package `ever-given` with `pip install ever-given`. The `ever-given` package mimics how we will run you container when uploaded to the website. For more information please see [`tutorials/ever_givenUsage.md`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/ever_givenUsage.md)


In your activated virtual environment, run the following:
* Batched: 
	```
	evergiven {container_name:tag} --container_type {type} --container_engine {engine} --file-receptor {path_to_receptor} --smiles {path_to_smiles_csv} --c_x {float} --c_y {float} --c_z {float} --sz_x {float} --sz_y {float} --sz_z {float}
	```
	```
	evergiven osatom/nanoluc-dock --container_type docker --container_engine docker --file-receptor {path_to_SAMPL-containers}/virtual_screening/data/5ibo.pdb --smiles {path_to_SAMPL-containers}/virtual_screening/data/smiles.csv --c_x 51.3155 --c_y 42.0975 --c_z 19.9691 --sz_x 40.307 --sz_y 38.3611 --sz_z 54.9179
	```
* Un-batched:
	```
	evergiven {container_name:tag} --container_type {type} --container_engine {engine} --file-receptor {path_to_receptor} --smiles {quoted_smiles_str} --c_x {float} --c_y {float} --c_z {float} --sz_x {float} --sz_y {float} --sz_z {float}
	```
	```
	evergiven osatom/nanoluc-dock --container_type docker --container_engine docker --file-receptor {path_to_SAMPL-containers}/virtual_screening/data/5ibo.pdb --smiles "ClC1=CC\2=C(OC3=C(C=CC=C3)\N=C2N4CCNCC4)C=C1" --c_x 51.3155 --c_y 42.0975 --c_z 19.9691 --sz_x 40.307 --sz_y 38.3611 --sz_z 54.9179
	```
* Please note:
	* `--container_type`: 
		* `docker` if container was built using `docker build` 
		* `singularity` if built container was built using `singularity build` or using a container sif file
	* `--container_engine`: the engine you will use to run your container `docker` or `singularity`
