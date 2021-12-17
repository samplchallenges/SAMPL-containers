# Conda Environment Instructions
In this document, we will explain how to create a conda environment inside your container. 

Containers are built/compiled using a set of instructions specified by a file called a [`Dockerfile`](https://docs.docker.com/engine/reference/builder/#:~:text=to%20Docker%20Hub.-,BuildKit,-%F0%9F%94%97). This instruction set can be shortened by starting from a "parent image". If you specify a parent image in your `Dockerfile`, your container will inherit everything setup and installed inside the parent image. 

In this document, we will show you two different ways to incorporate `miniconda` inside your container. The first section will inherit from a `miniconda` parent image, meaning we can use `conda` commands without any additional steps. The second section will describe what to do when the `miniconda` parent image is incompatible with your programs/workflow. In this case, we will need to explicitly specify the instructions to install miniconda in our `Dockerfile` build instructions. 

As a note, for SAMPL challenges you should always install your required python packages into the `base` environment. Please **DO NOT** create a separate virtual environment as this will cause command line parsing issues for our SMILES inputs 

## Section 1: Creating a Conda Environment Inside a Container using a `miniconda` Parent Image
### Part 1: Create your conda environment
> For this section, please ensure you have the Docker Desktop app and the Docker SDK installed. Please also ensure your Docker Desktop app is started up otherwise you will get a docker daemon error.
1. Start running a miniconda container
   * command: `docker run -it --rm continuumio/miniconda3`
2. Create a conda environment
   * command: `conda create --name <env_name> python=<version>`
3. Activate the conda environment
   * command: `conda activate <env_name>`
4. Install the packages your docking program will require
   * command: `conda install -c <channel> <package>`
   * command: `pip install <package>`
5. Test your environment to ensure all imports work correctly
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59) 
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import package
   >>> 
   >>> quit()
   ```
6. Once you are confident all packages have been installed correctly, export your environment to a yaml file. The output should look something like the code block below
   * command: `conda env export -n <env_name>`
   ```
   name: env_name
   channels:
     - conda-forge
     - defaults
   dependencies:
     - _libgcc_mutex=0.1=conda_forge
     ...
     - zstd=1.5.0=ha95c52a_0
   prefix: /opt/conda/envs/env_name
   ```
7. Copy the output from the export command in step 6 to be pasted into a file in step 10.
8. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
   * command `exit`
9. In the directory housing your container code, files, etc. create a file called `environment.yml` 
10. Open `environment.yml` for editing and paste in the output from step 6.
11. Change the first line of the file from `name: env_name` to `name: base`
12. Delete the last line of the file: `prefix: /opt/conda/envs/advenv`
    * The edits to `environment.yml` should now look similar to the first and last lines of the figure below:
      ![yaml](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/images/conda_env_yml.png)
14. Save and close `environment.yml`

### Part 2: Add the Build steps to install your environment
1. In the directory housing your container `environment.yml`, code, etc., create a file called `Dockerfile` (if you haven't already), and add the following lines to your `Dockerfile`. 
    * Please note we are using `conda env update` rather than `conda create` so your `base` environment is updated. 
    * Please **DO NOT** change this as using a conda environment other than `base` may prevent your program from properly accepting our command line input arguments.
    ```
    FROM continuumio/miniconda3:4.9.2-alpine  
    # tells the container to inherit from a miniconda container

    WORKDIR /opt/app/   
    # set the work directory

    COPY . ./    
    #  copy all the files and directories into the container

    RUN conda env update -f environment.yml && \
        conda clean --all --yes      
    # install the packages in environment.yml into containers
    ```
2. Save and close `Dockerfile` 
3. Build your your container

### Part 3: Test your container environment
> For this section, please ensure you have the [Docker Desktop](https://www.docker.com/products/docker-desktop) app and the [Docker SDK](https://pypi.org/project/docker/) installed. Please also ensure your Docker Desktop app is started up otherwise you will get a docker daemon error. 
1. In the directory housing your container `environment.yml`, `Dockerfile`, etc. build your container. Please ensure that if you have added an `ENTRYPOINT` line to your `Dockerfile` you comment out the `ENTRYPOINT` line (comments begin with `#` similar to Python). 
    * command: `docker build -t container-name:tag .`
2. Run your container
    * command: `docker run -it --rm container-name:tag`
3. Your command prompt should change after running this command. (See the code block below) If the command prompt doesn't change, sometimes the `ctrl-D` command remedies this.
    ```
    megosato ~ % docker run -it --rm container-name:tag
    root@7807d633195f:/opt/app#
    ```
4. Start up the Python interpreter
    * command: `python`
    ```
    root@7807d633195f:/opt/app# python
    Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59) 
    [GCC 7.5.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> 
    ```
 5. In the Python interpreter, `import` your required packages, to ensure they were installed properly and there are no errors.
    ```
    >>> import rdkit
    >>> 
    ```
6. Once you are confident all the packages have been properly installed, quit the python interpreter
    ```
    >>> quit()
    root@7807d633195f:/opt/app# 
    ```
7. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
    * command: `exit`



## Section 2: Creating a Conda Environment using Other Parent Images (besides `miniconda` image)
> In some cases, the `continuumio/miniconda` environment may not be compatible with the programs required by your docking protocol. In this case, you may need to use a different base image when building your container. Please see _______ to learn more about using other base images.
1. Go to [dockerhub](https://hub.docker.com/) and use the search bar to search for a container that meets your needs. Please save the name of the image you will use as your base, we will use it in later steps. 
	* For example, if I needed a container with a `gcc` compiler I would search for `gcc`, choose an image and locate the image name:
		![searchbar](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/images/dockerhub_search.png)
2. Once you have the name of the image you will use as your base, let's call it `image-to-use`, create a [Dockerfile](https://github.com/samplchallenges/SAMPL-containers/blob/3ddb358e64aa542c230da0af686d2fa3186108a9/tutorials/templates/docking/Dockerfile#L1) (if you have not already) and add the following: 
	* Please ensure that you replace `image-to-use` on the `FROM` line with your chosen base image from step 1. 
	```
	FROM image-to-use
	# tells the container to inherit from your chose container image

	WORKDIR /opt/app/   
	# set the work directory

	COPY . ./    
	#  copy all the files and directories into the container

	```
4. To install miniconda into your container at build time, we have found the following to the end of your `Dockerfile`:
	* Please note that we have found the following command sequence to be the most successful for installing miniconda into a container. However, this may not work for all parent images; some troubleshooting may be required if this is the case. 
	```
	RUN wget \
	    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	    && mkdir /root/.conda \
	    && bash Miniconda3-latest-Linux-x86_64.sh -b \
	    && rm -f Miniconda3-latest-Linux-x86_64.sh
	```

5. Once you have done the above steps, please follow the steps in [Section 1, Part 1](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/CondaEnvironmentInstructions.md#part-1-create-your-conda-environment) to finish create your `environment.yml` file. 
6. Upon creating the `environment.yml` file, add the following to the end of your `Dockerfile`
	* Please note we are using `conda env update` rather than `conda create` so your `base` environment is updated. 
	* Please **DO NOT** change this as using a conda environment other than `base` may prevent your program from properly accepting our command line input arguments.
	```
	RUN conda env update -f environment.yml && \
	    conda clean --all --yes      
	# install the packages in environment.yml into containers
	```
7. Finally, test your container conda environment using the steps in [Section 1, Part 3](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/CondaEnvironmentInstructions.md#part-3-test-your-container-environment)



## Section 3: What to do when Conda Environments are not Compatible
In some cases, miniconda3 may not be compatible with your programs inside your container or may be a source of licensing issues. If this is the case for you, this section details a work around we have found that only uses Python3 and `pip`.
1. Using the instructions in Section 2, step 1 and 2, you will need to find a container that fits your needs and has a `gcc` compiler, and then create a `Dockerfile`
2. In your `Dockerfile` from step 1, add the steps to install python. The following installs Python v3.6.0, but feel free adapt the following for a different Python version:
	```
	RUN apt-get update -y
	RUN apt-get install -y sudo build-essential
	
	RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
		tar xvf Python-3.6.0.tgz && \
		cd Python-3.6.0 && \
		./configure && \
		make -j 8 && \
		sudo make altinstall && \
		cd /opt/app/ && \
		rm -rf Python-3.6.0 && \
		rm Python-3.6.0.tgz
		# install python 3.6 for debian 8
		
	ENV LC_ALL C.UTF-8
	ENV LANG C.UTF-8
	# set environment variables required to use pip
	```
3. Add all the `pip` installable packages you require to the `install_requires` list in the [`setup.py`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/templates/docking/setup.py) file. 
4. Add the following to the end of your Dockerfile
	```
	RUN pip3.6 install .
	# pip install any extra python packages you will need
	```
5. Build your container to ensure there are no build issues
	* command: `docker build -t <container-name>:<tag> .`
6. If step 5 is successful, add the `ENTRYPOINT` for your program to your `Dockerfile`. In the line below, replace the `"entrypoint-from-setup.py"` with the entrypoint you set up in the `entry_points` string in the [`setup.py`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/templates/docking/setup.py) file. After this step, you should be good to finish developing your container
	```
	ENTRYPOINT ["entrypoint-from-setup.py"]
	```
7. If step 5 was unsuccessful, please continue with the following steps.
8. Remove the `RUN pip3.6 install .` line from your `Dockerfile`
9. Add the following line in place of `RUN pip3.6 install .`. Please add any other `pip` installable packages you need to the following line as well. 
	```
	RUN pip3.6 install click
	```
10. Archive or remove the `setup.py` file from your container directory as it is no longer needed for this step. 
11. Add the `ENTRYPOINT` build instruction below to your `Dockerfile`, you will need to replace `main.py` with the name of the file containing your containers main function. 
	```
	ENTRYPOINT ['python3.6', 'main.py']
	```
11. At the bottom of your the file containing your containers main function, you will need to add the following lines, replacing `main_function()` with your main function. You do not need to pass any arguments to the `main_function` as these will be handled by `click`.
	```
	if __name__ == "__main__":
		main_function()
	```
12. Build your container to ensure there are no build issues.
	* command: `docker build -t <container-name>:<tag> .`
	
13. If your container builds correctly, run your container with the `--help` flag. The output should look similar to the following code block
	* command: `docker run -it --rm <container-mane>:<tag> --help`
	```
	(base) megosato@Admins-MacBook-Pro tutorial % docker run -it --rm osatom/adv-tutorial --help
	Usage: main.py [OPTIONS]

	Options:
	  --receptor PATH      path of receptor PDB to dock the ligand into
			       [required]
	  --smiles TEXT        SMILES str of ligand to be docked. quote to prevent CLI
			       errors "CCC"  [required]
	  --hint PATH          path of hint ligand complex for docking region hint
			       [required]
	  --hint-molinfo TEXT  residue name of the ligand in the hint complex
			       [required]
	  --hint-radius FLOAT  box size of the box to dock into  [required]
	  --output-dir TEXT    Output directory for receptor and docked_ligand files
	  --debug              prints debug print statements when --debug flag is used
	  --help               Show this message and exit.
	 ```
