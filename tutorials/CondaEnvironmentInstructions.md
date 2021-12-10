# Conda Environment Instructions
WHAT IS A BASE IMAGE


## Section 1: Creating a Conda Environment Inside a Container using a `miniconda` Base Image
### Part 1: Create your conda environment
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
    * The edits to `environment.yml` should now look similar to the figure below:
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



## Creating a Conda Environment using Other Base Images (besides `miniconda` image)
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
4. To install miniconda into your container at build time, add the following to your `Dockerfile`:
	```
	RUN wget \
	    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	    && mkdir /root/.conda \
	    && bash Miniconda3-latest-Linux-x86_64.sh -b \
	    && rm -f Miniconda3-latest-Linux-x86_64.sh
	```

5. Once you have done the above steps, please follow the steps in [Section 1, Part 1] to finish create your `environment.yml` file. 
6. Upon creating the `environment.yml` file, add the following to the end of your `Dockerfile`
	* Please note we are using `conda env update` rather than `conda create` so your `base` environment is updated. 
	* Please **DO NOT** change this as using a conda environment other than `base` may prevent your program from properly accepting our command line input arguments.
	```
	RUN conda env update -f environment.yml && \
	    conda clean --all --yes      
	# install the packages in environment.yml into containers
	```
7. Finally, test your container conda environment using the steps in [Section 1, Part 3]


* Any additional steps to install your required programs should also be added to your `Dockerfile` for more information, please see [Section 1.4](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#14-download-and-prepare-the-command-line-programs-autodock-vina-and-mgl-tools-executables-for-use-in-the-docking-container) and [Section 1.5](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#15-install-autodock-vina-and-mgl-tools-into-your-container) of the tutorial. 

