# Conda Environment Instructions
In this document, we will explain how to create a conda environment inside your container using both Docker and Singularity.

Containers are built/compiled using a set of instructions specified by a file . This instruction set can be shortened by starting from a "parent image". If you specify a parent image, your container will inherit everything setup and installed inside the parent image.

In this document, we will show you two different ways to incorporate `miniconda` inside your container. The first section will inherit from a `miniconda` parent image, meaning we can use `conda` commands without any additional steps. The second section will describe what to do when the `miniconda` parent image is incompatible with your programs/workflow. In this case, we will need to explicitly specify the instructions to install miniconda in our instruction file build instructions.

As a note, for SAMPL challenges you should always install your required python packages into the `base` environment. Please **DO NOT** create a separate conda virtual environment as this will cause command line parsing issues for our SMILES inputs.

called a [`Dockerfile`](https://docs.docker.com/engine/reference/builder/#:~:text=to%20Docker%20Hub.-,BuildKit,-%F0%9F%94%97)

# Docker Engine Conda Environment Creation
## Section 1: Creating a Conda Environment Inside a Container using a `miniconda` Parent Image
### Part 1: Create your conda environment
> For this section, please ensure you have the Docker Desktop app and the Docker SDK installed. Please also ensure your Docker Desktop app is started up otherwise you will get a docker daemon error.

1. Start running an interactive miniconda container instance
   * command: `docker run -it --rm continuumio/miniconda3`
2. Install the python version necessary
   * command: `conda install python=<version>`
3. Install the packages your docking program will require into the base environment
   * command: `conda install -c <channel> <package>`
   * command: `pip install <package>`
4. Test your environment to ensure all imports work correctly
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import <module>
   >>>
   >>> quit()
   ```
5. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
   * command: `exit`

### Part 2: Add the Build steps to install your environment
1. In the directory housing your container code, etc., create a file called `Dockerfile` (if you haven't already), and add the following lines to your `Dockerfile`.
    * Please note we are using `conda env update` rather than `conda create` so your `base` environment is updated.
    * Please **DO NOT** change this as using a conda environment other than `base` may prevent your program from properly accepting our command line input arguments.
    ```
    FROM continuumio/miniconda3:4.9.2-alpine  
    # tells the container to inherit from a miniconda container

    WORKDIR /opt/app/   
    # set the work directory

    COPY . ./    
    #  copy all the files and directories into the container


    RUN conda update conda && \
        conda install python=<version> && \
        conda install -c <channel> <package> && \
        ...
        conda clean --all --yes
    # install required conda packages      

    RUN pip install <package>
    # install required pip packages
    ```
2. Save and close `Dockerfile`

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
    >>> import <module>
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
5. Build your container.
  * command: `docker build test-env .`
6. Start running an interactive test_env container instance
  * command: `docker run -it --rm test-env`
7. Install the python version necessary
   * command: `conda install python=<version>`
8. Install the packages your docking program will require into the base environment
   * command: `conda install -c <channel> <package>`
   * command: `pip install <package>`
9. Test your environment to ensure all imports work correctly
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import <module>
   >>>
   >>> quit()
   ```
10. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
   * command `exit`
11. To install your environment into your container, follow the steps in [Section 2, Part  2]()
12. Finally, test your container conda environment using the steps in [Section 1, Part 3](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/CondaEnvironmentInstructions.md#part-3-test-your-container-environment)

## Section 3: What to do when Conda Environments are not Compatible
> In some cases, miniconda3 may not be compatible with your programs inside your container or may be a source of licensing issues. If this is the case for you, this section details a work around we have found that only uses Python3 and `pip`.

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
3. Add all the `pip` installable packages you require to the `RUN` commands in the Dockerfile
	```
	RUN pip3.6 install package1 package2 ...
	# pip install any extra python packages you will need
	```
4. Build your container to ensure there are no build issues
	* command: `docker build -t <container-name>:<tag> .`
5. If step 4 is successful, add the `ENTRYPOINT` for your program to your `Dockerfile`. In the line below, replace the `"entrypoint-from-setup.py"` with the entrypoint you set up in the `entry_points` string in the [`setup.py`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/templates/docking/setup.py) file. After this step, you should be good to finish developing your container
	```
	ENTRYPOINT ["entrypoint-from-setup.py"]
	```
6. If step 5 was unsuccessful, please continue with the steps 7 to 10.
7. Archive or remove the `setup.py` file from your container directory as it is no longer needed for this step.
8. Add the `ENTRYPOINT` build instruction below to your `Dockerfile`, you will need to replace `main.py` with the name of the file containing your containers main function.
	```
	ENTRYPOINT ['python3.6', 'main.py']
	```
9. At the bottom of your the file containing your container's main function, you will need to add the following lines, replacing `main_function()` with your main function. You do not need to pass any arguments to the `main_function` as these will be handled by `click`.
	```
	if __name__ == "__main__":
		main_function()
	```
10. Build your container to ensure there are no build issues.
	* command: `docker build -t <container-name>:<tag> .`

11. If your container builds correctly, run your container with the `--help` flag. The output should look similar to the following code block
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

# Singularity Engine Conda Environment Creation
## Section 1: Creating a Conda Environment Inside a Container using a `miniconda` Parent Image
### Part 1: Create your conda environment
> For this section, please ensure you have [Singularity 3.7.2 or 3.7.3](https://sylabs.io/guides/3.7/user-guide/quick_start.html).

1. Start running an interactive miniconda container instance
   * command: `singularity run docker://continuumio/miniconda3`
2. Install the python version necessary
   * command: `conda install python=<version>`
3. Install the packages your docking program will require into the base environment
   * command: `conda install -c <channel> <package>`
   * command: `pip install <package>`
4. Test your environment to ensure all imports work correctly
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import <module>
   >>>
   >>> quit()
   ```
5. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
   * command `exit`

### Part 2: Add the Build steps to install your environment
1. In the directory housing your container code, etc., create a file called `buildfile` (if you haven't already), and add the following lines to your `buildfile`.
    * Please note we are using `conda env update` rather than `conda create` so your `base` environment is updated.
    * Please **DO NOT** change this as using a conda environment other than `base` may prevent your program from properly accepting our command line input arguments.
    ```
    Bootstrap: docker
    From: continuumio/miniconda3

    %files

    %post
    conda update conda && \
        conda install python=<version> && \
        conda install -c <channel> <package> && \
        ...
        conda clean --all --yes
    # install required conda packages      

    pip install <package>
    # install required pip packages
    ```
2. Save and close `buildfile`

### Part 3: Test your container environment
> For this section, please ensure you have the [Docker Desktop](https://www.docker.com/products/docker-desktop) app and the [Docker SDK](https://pypi.org/project/docker/) installed. Please also ensure your Docker Desktop app is started up otherwise you will get a docker daemon error.

1. In the directory housing your container `environment.yml`, `Dockerfile`, etc. build your container. Please ensure that if you have added an `ENTRYPOINT` line to your `Dockerfile` you comment out the `ENTRYPOINT` line (comments begin with `#` similar to Python).
    * command: `singularity build --fakeroot container_name.sif buildfile`
2. Run your container
    * command: `singularity run container_name.sif`
3. Your command prompt should change after running this command. (See the code block below) If the command prompt doesn't change, sometimes the `ctrl-D` command remedies this.
    ```
    vagrant@ubuntu-bionic:adv-tut-singularity-base$ singularity shell
    INFO:    Using cached SIF image
    Singularity>
    ```
4. Start up the Python interpreter
    * command: `python`
    ```
    Singularity> python
    Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
    [GCC 7.5.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>
    ```
5. In the Python interpreter, `import` your required packages, to ensure they were installed properly and there are no errors.
    ```
    >>> import <module>
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
  Bootstrap: docker
  From: image-to-use
	# tells the container to inherit from your chose container image

  %files

  %post
  conda update conda && \
      conda install python=<version> && \
      conda install -c <channel> <package> && \
      ...
      conda clean --all --yes
  # install required conda packages      

  pip install <package>
  # install required pip packages
	```
4. To install miniconda into your container at build time, we have found the following to the `%post` section of your `buildfile`:
	* Please note that we have found the following command sequence to be the most successful for installing miniconda into a container. However, this may not work for all parent images; some troubleshooting may be required if this is the case.
	```
  %post
	wget \
	    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
	    && mkdir /root/.conda \
	    && bash Miniconda3-latest-Linux-x86_64.sh -b \
	    && rm -f Miniconda3-latest-Linux-x86_64.sh
	```
5. Build your container.
  * command: `singularity build --fakeroot test-env buildfile`
6. Start running an interactive test_env container instance
  * command: `singularity shell test-env`
7. Install the python version necessary
   * command: `conda install python=<version>`
8. Install the packages your docking program will require into the base environment
   * command: `conda install -c <channel> <package>`
   * command: `pip install <package>`
9. Test your environment to ensure all imports work correctly
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import <module>
   >>>
   >>> quit()
   ```
10. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
   * command `exit`
11. To install your environment into your container, follow the steps in [Singularity Section 1, Part  2]()
12. Finally, test your container conda environment using the steps in [Singularity Section 1, Part 3]()

## Section 3: What to do when Conda Environments are not Compatible
> In some cases, miniconda3 may not be compatible with your programs inside your container or may be a source of licensing issues. If this is the case for you, this section details a work around we have found that only uses Python3 and `pip`.

1. Using the instructions in Section 2, step 1 and 2, you will need to find a container that fits your needs and has a `gcc` compiler, and then create a `buildfile`
2. In your `buildfile` from step 1, add the steps to install python. The following installs Python v3.6.0, but feel free adapt the following for a different Python version:
	```
  %post
  apt-get update -y
	apt-get install -y sudo build-essential

	wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
		tar xvf Python-3.6.0.tgz && \
		cd Python-3.6.0 && \
		./configure && \
		make -j 8 && \
		sudo make altinstall && \
		cd /opt/app/ && \
		rm -rf Python-3.6.0 && \
		rm Python-3.6.0.tgz
		# install python 3.6 for debian 8

  LC_ALL=C.UTF-8
	LANG=C.UTF-8
	# set environment variables required to use pip
	```
3. Add all the `pip` installable packages you require to the `%post` commands in the `buildfile`
	```
	pip3.6 install package1 package2 ...
	# pip install any extra python packages you will need
	```
4. Build your container to ensure there are no build issues
	* command: `singularity build --fakeroot container-name.sif buildfile`
5. If step 4 is successful, add the `%runscript` for your program to your `buildfile`. In the line below, replace the `"entrypoint-from-setup.py"` with the entrypoint you set up in the `entry_points` string in the [`setup.py`](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/templates/docking/setup.py) file. After this step, you should be able to finish developing your container. Please skip to step 11.
	```
  %runscript
  exec <entrypoint-from-setup.py> $@
	```
6. If step 5 was unsuccessful, please continue with the steps 7 to 10.
7. Archive or remove the `setup.py` file from your container directory as it is no longer needed for this step.
8. Replace the `%runscript` build instruction below to your `buildfile`, you will need to replace `main.py` with the name of the file containing your containers main function.
	```
	exec python3.6 main.py $@
	```
9. At the bottom of your the file containing your container's main function, you will need to add the following lines, replacing `main_function()` with your main function. You do not need to pass any arguments to the `main_function` as these will be handled by `click`.
	```
	if __name__ == "__main__":
		main_function()
	```
10. Build your container to ensure there are no build issues.
	* command: `singularity build --fakeroot container-name.sif buildfile`

11. If your container builds correctly, run your container with the `--help` flag. The output should look similar to the following code block
	* command: `singularity run container-name.sif --help`
	```
	(base) megosato@Admins-MacBook-Pro tutorial % singularity run container-name.sif --help
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
