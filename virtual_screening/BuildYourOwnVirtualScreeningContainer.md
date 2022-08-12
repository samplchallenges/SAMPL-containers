# Building Your Own Virtual Screening Container
> This document details requirements and tips for writing a docking container for SAMPL-challenges. For Python template files that follow this guide, please see [SAMPL-containers/virtual_screening/templates/docker](https://github.com/samplchallenges/SAMPL-containers/tree/main/virtual_screening/templates/docker). For an example of a run-able virtual screening main file, please see [SAMPL-containers/virtual_screening/examples/adv-screen-docker/main.py](https://github.com/samplchallenges/SAMPL-containers/blob/main/virtual_screening/examples/adv-screen-docker/main.py). This guide is written under the assumption the reader has already gone through the [Docking Tutorial](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md). 
> 
> 
> **IMPORTANT: Your container should be failsafe. If your program(s) fail on a specific output, you should still present a result/prediction. Perhaps you can predict a default value of 0 (non-binder) for whether the compound binds, and a value of 0 for the score. For ideas on how to make your container failsafe, please see [WritingFailsafeContainers.md](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/WritingFailsafeContainers.md).**

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


## Program Logs
* Any output to `stdout` or `stderr` will be logged with timestamps associated with each output. These logs will be made accessible to you.
* Please print general logging info to `stdout` and error messages to `stderr` as is convention.
* As stated in [OutputRequirements](https://github.com/samplchallenges/SAMPL-containers/blob/main/virtual_screening/BuildYourOwnVirtualScreeningContainer.md#output-requirements), the last two lines of `stdout` output must be your two `key value` pairs. 
* Please remember to [flush your print statements](https://www.geeksforgeeks.org/python-output-using-print-function/) so your output will be appear as it is executed


## A Note about Intermediate Files
* Please do not use `output_dir` to store your intermediate files
* You can store your intermediate files in a temporary directory or in any directory other than the `output_dir`
* Because we will be running your containers with batches of inputs, rather than 1 input at a time, it may be good to delete files you no longer need as you go to avoid storage/memory issues. 


## Including your own Python Modules
If you modularize your code and include your own python modules, you will need to follow the steps below. For an example with using extra python modules beyond just main.py, please see [SAMPL-containers/virtual_screening/examples/adv-tutorial](https://github.com/samplchallenges/SAMPL-containers/tree/main/pose_prediction/examples/adv-docker).
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
		```
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
    * Docker: [`SAMPL-league/virtual_screening/examples/adv-base-docker/Dockerfile`](https://github.com/samplchallenges/SAMPL-containers/tree/main/virtual_screening/examples/adv-base-docker/Dockerfile)
    * Singularity: WIP
* To run a command line program from within a Python module, consider using the [`subprocess`](https://docs.python.org/3/library/subprocess.html) library or similar from the Python3 library
    * Please see [`SAMPL-league/virtual_screening/examples/adv-screen-docker/autodock.py`](https://github.com/samplchallenges/SAMPL-containers/tree/main/virtual_screening/examples/adv-screen-docker/autodock.py)
