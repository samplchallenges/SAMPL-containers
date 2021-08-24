# Introduction to SAMPL Containerized Methods

## Purpose:
In [SAMPL4](https://link.springer.com/article/10.1007%2Fs10822-013-9702-2), we learned that human knowledge can be a key factor influencing the success of a computational drug discovery method. To work around this finding, we are creating an automated arm of SAMPL challenges to run methods head-to-head without human intervention. To accomplish this, we will use containerized [Docker](https://www.docker.com/resources/what-container) methods.

The following tutorial is meant to teach the basics of building a simple container to predict docking poses. Here, we use both Python code and command-line programs (specifically, Autodock Vina and MGL Tools).

## Background and Prerequisites

### Important Note on "Docker" versus "Docking":
   Please note that "Docker" and "docking" are two separate things.
* **"Docker"** is a program that allows you to containerize methods, essentially allowing you to distribute your method or approach in a reusable way where it can be used reproducibly without human intervention
* **"Docking"** describes predicting the structure of a complex, in this case a protein-ligand complex.

This terminology is unfortunately not something we can change.

### Expected Background Knowledge
* Basic knowledge of Python
* Basic knowledge of Linux/UNIX command line


### Software Requirements
* Linux or Mac operating system
* [Docker Desktop](https://www.docker.com/products/docker-desktop)
* [Docker SDK for Python](https://pypi.org/project/docker/)
* [Python 3](https://www.python.org/downloads/)

## Getting started

### Brief Docker Usage Tips
* To build an image, ensure you are in the directory with your Dockerfile and container code, then run `docker build -t <name>:<tag/version> .`
  * Examples: `docker build -t adv:0.1 .` or `docker build -t adv:latest .`
* Use the command `docker images` to list out your built images
* To delete Docker images, use `docker images` to list your current images and their IMAGE IDs, then run the command `docker image rm <IMAGE IDs>`


### Pre-Built Autodock Vina Container

A working version of the Autodock Vina container we will build in this tutorial can be found at [Docker Hub under osatom/adv-tutorial](https://hub.docker.com/repository/docker/osatom/adv-tutorial). To play with this container, please use the following steps:
1. Pull the "adv-tutorial" docker container: `docker pull osatom/adv-tutorial:latest`
2. Change directories into the "examples" directory: `cd examples`
3. Run the command: `python ever_given/run.py osatom/adv-tutorial:latest --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
4. The results will be stored in the directory "examples/evergiven_output"


# Tutorial: Build an AutoDock Vina Containerized Method

## Outline:
* [Section 1: Build the Autodock Vina base container](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#section-1-build-the-autodock-vina-base-container)
   * [1.1 Setup](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#11-setup)
   * [1.2 Starting a docker file and creating a conda environment](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#12-create-a-conda-environment)
   * [1.3 Install conda environment (from Section 1.2) into your container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#13-create-a-dockerfile)
   * [1.4 Add the Autodock Vina and MGL Tools executables](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#14-add-the-autodock-vina-and-mgl-tools-executables)
   * [1.5 Update the Dockerfile to include Autodock Vina and MGL Tools installation steps](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#15-update-the-dockerfile-to-include-autodock-vina-and-mgl-tools-installations)
   * [1.6 Build the base container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#16-build-the-base-container)
* [Section 2: Build the Autodock Vina Docing methods container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#section-2-build-the-container-with-autodock-vina-docking-methods)
   * [2.1 Setup](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#21-setup)
   * [2.2 Add the docking code](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#22-add-the-docking-code)
   * [2.3 Create a setup.py file](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#23-create-a-setuppy-file)
   * [2.4 Create a Dockerfile](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#24-create-a-dockerfile)
   * [2.5 Build the docking container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#24-create-a-dockerfile)
 * [Section 3: Test/Run your container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#section-3-testrun-your-container)

## Section 1: Build the Autodock Vina base container
In this section, we will build a base container that has all necessary packages and programs installed. This way, as we write our docking code it will build quickly since only the docking code will need to be built.


### 1.1: Setup
1. Create a directory called "adv-tutorial-base" in the examples directory
   * command: `mkdir adv-tutorial-base`
2. Change directories to "adv-tutorial-base"
   * command: `cd adv-tutorial-base`


### 1.2: Starting a pre-made docker container and creating a conda environment

> In 1.2, we will start the pre-made minconda container, "continuumio/miniconda3", which contains a pre-installed version of miniconda. This will allow us to interactively create the conda environment we will use inside the container. Even if you have a conda environment installed locally, you will need to complete this step. The container is an isolated virtual machine, so it will not have access to your local conda environment. 
> 
> Because our container will inherit from the miniconda container, any environment we create while interatively using the miniconda container should install without additional issues. Building a conda environment outside the miniconda container often results in multiple rounds of trial and error and incompatible packages. For more detailed/generalized instructions please see [CondaEnvInstructions.pdf](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/CondaEnvInstructions.pdf)
1. Start up Docker Desktop, which will start the Docker daemon. If this is your first time starting Docker Desktop, the application may need to complete some additional installation steps or updates. 
2. Open a terminal
3. Start the container. Upon running this command your command line prompt should change. This means you are now inside the container interacting with it command line and contents. 
   * command: `docker run -it --rm continuumio/miniconda3`
   ```
   megosato@Admins-MacBook-Pro ~ % docker run -it --rm continuumio/miniconda3
   root@7f02be71557e:/# 
   ```
4. Create a conda env called "advenv"
   * command: `conda create --name advenv`
5. Activate advenv: `conda activate advenv`
6. Install rdkit
   * command: `conda install -c conda-forge rdkit`
7. Install mdtraj
   * command: `conda install -c conda-forge mdtraj`
8. Export the environment. Upon running this command, an organized list of the environment packages will be printed out to your console (something like the code block below).
   * command: `conda env export -n advenv`
   ```
   (advenv) root@7f02be71557e:/# conda env export -n advenv
   name: advenv
   channels:
     - conda-forge
     - defaults
   dependencies:
     - _libgcc_mutex=0.1=conda_forge
     ...
     - zstd=1.5.0=ha95c52a_0
   prefix: /opt/conda/envs/advenv
   ```
9. Copy the output from the export command in step 7 to be pasted into a file in step 11.
10. Exit the container
   * command: `exit`
11. Create and open a file called "environment.yml" and paste the output you previously copied at step 8
12. Change the first line of the file `name: advenv` to `name: base`
13. Delete the last line of the file: `prefix: /opt/conda/envs/advenv`
14. Save the changes to environment.yml and exit


### 1.3: Install conda environment (from [Section 1.2](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#12-starting-a-pre-made-docker-container-and-creating-a-conda-environment) into your container  

> In 1.3, we will begin creating a Dockerfile which contains the instructions required to build the base container. For more information on how to visit a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image).
1. Create and open a file called "Dockerfile"
2. Copy the following lines into Dockerfile
   ```
   FROM continuumio/miniconda3:4.9.2-alpine  
   # tells the container to inherit from a miniconda container

   WORKDIR /opt/app/   # set the work directory

   COPY . ./    #  copy all the files and directories into the container

   RUN conda env update -f environment.yml && \
       conda clean --all --yes      # install the packages in environment.yml into containers

   ENV PATH="/root/.local/bin:$PATH"      # set the path
   ```
3. Save the changes to Dockerfile and exit
4. Build your container to ensure there are no build issues, so far. 
   * command: `docker build -t adv-tutorial-base-test .`
5. If your build from the previous step (step 4) completed without issue, please move on to Section 1.4, otherwise some troubleshooting may be necessary.


### 1.4: Add the Autodock Vina and MGL Tools executables

> In 1.4, we will incorporate the command line tools Autodock Vina and MGL Tools into our base container. Please do not change which installers you download based on your native operating system (OS) because the OS used inside docker container, Linux x86, may be different than your native OS. For example, I am currently on a Mac, but the OS inside the docker container is Linux x86, so any installers that work for my native Mac OS would not work inside my docker container. 
> 
> When building your own container, this is where you would add in any command line program files.
1. Create a directory called dependencies
   * command: `mkdir dependencies`
2. Download Autodock Tools linux x86 "autodock_vina_1_1_2_linux_x86.tgz" from http://vina.scripps.edu/download.html
3. Untar "autodock_vina_1_1_2_linux_x86.tgz"
   * command: `tar -xvf dependencies/autodock_vina_1_1_2_linux_x86.tgz`
4. Delete the .tgz file:
   * command: `rm autodock_vina_1_1_2_linux_x86.tgz`
5. Rename "autodock_vina_1_1_2_linux_x86/" to "adv/"
   * command: `mv autodock_vina_1_1_2_linux_x86 adv`
6. Move the "adv" directory inside the "dependencies" directory
   * command `mv adv dependencies`
7. Download MGL Tools linux x86 `mgltools_x86_64Linux2_1.5.6.tar.gz` from http://mgltools.scripps.edu/downloads
8. Untar "mgltools_x86_64Linux2_1.5.6.tar.gz"
   * command `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
9. Delete the .tgz file:
   * command: `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
10. Rename "mgltools_x86_64Linux2_1.5.6.tar.gz" to "mgl"
      * command: `mv mgltools_x86_64Linux2_1.5.6.tar.gz mgl`
11. Move "mgl" directory inside the "dependencies" directory:
      * command: `mv mgl dependencies/`
12. Open "dependencies/mgl/install.sh"
13. Change line 6 from `TarDir=""` to `TarDir="/opt/app/dependencies/mgl/"`
14. Change line 7 from `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
15. Close and save "dependencies/mgl/install.sh"


### 1.5: Update the Dockerfile to include Autodock Vina and MGL Tools installation steps
> In 1.6, we will add in the installation commands for Autodock Vina and MGL Tools to the Dockerfile. When building your own container, this is where you would add in any command line program installation steps.
1. Open the Dockerfile and paste the following lines
   ```
   RUN /opt/app/dependencies/mgl/install.sh

   ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"

   RUN /opt/app/dependencies/adv/bin/vina --help
   ```
2. Save the changes to the `Dockerfile` and exit


### 1.6: Build the base container

> In 1.6, we will build a docker image using "docker build", which we will inherit from in the Section 2.

1. Build the base container
   * command: `docker build -t adv-tutorial-base .`


## Section 2: Build the container with Autodock Vina Docking methods

> In this section, we will build our docking container with runnable docking code.

### 2.1: Setup

1. Create a directory called "adv-tutorial" in the examples directory
   * command: `mkdir adv-tutorial`
2. Change directories to "adv-tutorial"
   * command: `cd adv-tutorial`

### 2.2: Add the docking code

> In 2.2, we will incorporate the docking code into our container directory. When building your own container, this is where you would add in your methods. However, for the sake of simplicity, we will be using pre-written docking code. Please see [examples/ContainerRequirements.md](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/ContainerRequirements.md) for more information on the inputs and kwargs required of each main function.
1. Copy the AutoDock class file from "examples/adv/autodock.py" to "adv-tutorial"
   * command: `cp ../adv/autodock.py .`
2. Copy the AutoDock main file from "examples/adv/main.py" to "adv-tutorial"
   * command: `cp ../adv/main.py .`

### 2.3: Create a setup.py file

> In 2.3, we will customize a setup.py file to to match the python modules we have written.
1. Create and open a file called "setup.py"
2. Copy and paste the following into setup.py
    ```
    from setuptools import setup

    setup(
        name='AutoDock-rdkit',
        version='0.1',
        py_modules=[
            ,
    	  ],
        install_requires=[
            ,
        ],
        entry_points='''
            [console_scripts]
        '''
    )
    ```
2. Modify the "py_modules" list, by adding the two modules with our docking code from the previous subsection ([2.2](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#22-add-the-docking-code)): "autodock" and "run_autodock"
    ```
    py_modules=[
        autodock,
        run_autodock,
    ]
    ```
3. Modify the "install_requires" list, adding "Click" a pip installable package we did not add to the previous base build
    ```
    install_requires=[
        Click,
    ]
    ```
4. Modify the "entry_points" section, by adding an entry point in the format `{command-to-call-in-Dockerfile}={py_module_with_main}:{function_to_run}`. The python module "main.py" contains the main function "main_function".
    ```
    entry_points='''
        [console_scripts]
        run-autodock=main:main_function
    '''
    ```
5. Save and close setup.py

### 2.4: Create a Dockerfile

> In 2.4, we will create a Dockerfile which contains the instructions required to build the docking container, as well as the entry_point (see [2.3](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#23-create-a-setuppy-file)) which tells the container which file and function to run. For more information on how to visit a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image).
1. Create and open a file called "Dockerfile"
2. Copy and paste the following into Dockerfile
    ```
    FROM

    WORKDIR /opt/app/

    COPY

    RUN pip install .

    ENV PATH="/root/.local/bin:$PATH"

    RUN ls -l /opt/app

    ENTRYPOINT []
    ```

2. Next to "FROM" add the name of the base build from Section 1
   * `FROM adv-tutorial-base`
3. Next to "COPY" add the names of all files necessary to run our docking program, including "setup.py", "autodock.py", and "main.py" as well as the container directory to copy them into "/opt/app" or "./"
   * `COPY setup.py autodock.py main.py ./`
4. Next to "ENTRYPOINT", add the "entry_point" you declared in step 4 of the previous subsection ([2.3](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#23-create-a-setuppy-file)) inside the brackets, in quotations.
   * `ENTRYPOINT ["run-autodock"]`


### 2.5: Build the docking container

> In 2.5, we will build a docker image that will execute our docking program when run. For more information about the `docker build`, please see the [official Docker documentation](https://docs.docker.com/engine/reference/commandline/build/).
1. Build the container
   * command: `docker build -t adv-tutorial .`


## Section 3: Test/Run your container
In this section, we will use the wrapper ever_given to run the docking container. ever_given mimics the infrastructure we will use to run your container on the SAMPL-league website, making it a great way to test that you container will run properly ahead of uploading to the website.
1. Change directories into "examples":
   * command: `cd examples`
2. Run the container
   * command: `python ever_given/run.py adv-tutorial --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
3. Examine your container outputs in the directory "examples/evergiven_output"


# Other Important Information
* For more detailed tutorials on how to use docker please see the following resources:
  * [Official Docker Documentation](https://docs.docker.com/get-started/)
  * [Brief Docker Tutorial (12m)](https://www.youtube.com/watch?v=YFl2mCHdv24)
  * [Docker Beginners Course (2hrs)](https://www.youtube.com/watch?v=fqMOX6JJhGo)
  * [Brief video on how Docker containers work (15m)](https://www.youtube.com/watch?v=rOTqprHv1YE)
  * YouTube is a great resource for learning Docker, feel free to search for other tutorials that suit your specific needs as well
* To upload your own container to docker.io create a [Docker Hub account](https://hub.docker.com/signup)
  * Please note that free Docker Hub accounts only receive 1 private repository.
  * If you have sensitive information inside your containers, please ensure you are using a private repository.
* Eventually, upon building enough Docker images, you may begin to run out of memory. Please remember to regularly delete any Docker images you no longer need. (See [Brief Docker Usage Tips](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#brief-docker-usage-tips): Bullet 3)
