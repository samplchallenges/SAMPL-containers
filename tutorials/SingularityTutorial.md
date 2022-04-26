# Introduction SAMPL Containerized Methods using the Singularity Container Engine

## Purpose
The following tutorial is meant to teach the basics of building a simple container to predict docking poses using the [Singularity](https://sylabs.io/guides/3.5/user-guide/introduction.html#:~:text=Why%20use%20containers%3F%C2%B6). Here, we use both Python code and command-line programs (specifically, Autodock Vina and MGL Tools).

## An Important Disclaimer
If you are developing your container/methods on a High Performace Computing (HPC) Cluster, you will more than likely need to use the Singularity container engine. Most HPC clusters will not have the Docker program installed. Please see this tutorial on how to build a container using Singularity.

If you have the ability to use either [Docker](https://www.docker.com/resources/what-container/) or Singularity as your container engine, please use Docker. The accepted best practice for containers is to build and store your containers as Docker images. A tutorial to build a Docker containerized method can be found [here](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md)

## Background and Prerequisites

### Expected Background Knowledge
* Basic knowledge of [Python3](https://www.python.org/download/releases/3.0/)
* Basic knowledge of [Linux/UNIX command line](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)

### Software Requirements
* Linux operating system or Linux Virtual Machine
* [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html)
* git command line tool


## Getting Started

### Brief Singularity Background

### Pre-Built Autodock Vina Singularity Container
A working version of the Autodock Vina container we will build in this tutorial with the Singularity engine is available in our google drive as the file [adv-tutorial.sif](https://drive.google.com/file/d/1yEKSLU6AKZzECUiTRpOoFVk6u-Bn44aY/view?usp=sharing). 
1. In a new virtual environment with python 3.8 install the ever-given package. Ensure the package version is the latest version by checking on [ever-given's PyPi](https://pypi.org/project/ever-given/)
      * command: `pip install ever-given`
2. If you haven't already, clone this repository
      * command: `git clone https://github.com/samplchallenges/SAMPL-containers.git`
4. Download the [adv-tutorial.sif](https://drive.google.com/file/d/1yEKSLU6AKZzECUiTRpOoFVk6u-Bn44aY/view?usp=sharing) using the link provided. 
5. Change directories into the "SAMPL-containers/tutorials" directory
      * command: `cd SAMPL-containers/tutorials`
6. Run the container
      * command: 
         ```
         evergiven {path_to_file}/adv-tutorial.sif \
            --container-engine singularity --container-type singularity_local \
            --file-receptor data/receptor.pdb \
            --file-hint data/hint.pdb --hint-radius 6 --hint-molinfo "E4Y" \
            --smiles "c1ccc(C(C)C)cc1CNCC(O)(O)[C@@H](NC(=O)[C@@H]2C)C[C@H](C)CCCCCCCCC(=O)N2C" \
            --output-keys docked_ligand,receptor
         ```
7. The results files will be stored in the directory `tutorials/evergiven_output`
      * Expected files:
           * prepped receptor file: `rec-dock.pdb`
           * docked ligand file: `best_dock.pdb`

# Tutorial: Build an AutoDock Vina Containerized Method using the Singularity Container Engine
> This tutorial is separated into six parts: (1) the virtual environment and dependency container build, (2) the docking container build, (3) docking using the Docker container, (4) trouble shooting, (5) building your own docking container, and (6) miscellaneous important information.
> 
> We have separated the base virtual environment build (conda environment and command-line programs) from the implementation build (docking) to illustrate container inheritance and to improve the time required to build the container. In Section 1, we will build an adv-tutorial-base container that includes all the environment and software installations necessary to run AutoDock Vina from a Python script. In Section 2, we will inherit the environment from the Section 1 adv-tutorial-base container to write and build our run-able docking code in the adv-tutorial container.
> 
> This inheritance scheme is also a good practice to improve the time required to build the docking container. When you are writing your docking container, most of the time, you will only make changes to the new docking code you have written and thus only need to re-build the new docking code, rather than re-build the new code AND the environment/software installations. By separating the docking container build from the environment container build, we can build just the docking container each time we change the docking code, improving build time. For example, building only new code (adv-tutorial) typically takes under 5 seconds, but building new code plus environment and software installations (adv-tutorial plus adv-tutorial-base) takes upwards of 115 seconds. This can be likened to independently compiling parts of a program to improve compile time as a whole later.


## Outline

## Section 1: Build the Autodock Vina base container
### 1.1: Setup
1. Open a terminal
2. If you have not already, in a new virtual environment with python 3.8 install the ever-given package. Please ensure the version you have installed the latest version listed on [ever-given's PyPi](https://pypi.org/project/ever-given/).
    * command: `pip install ever-given`
3. If you have not already, clone the SAMPL-containers repository
    * command: `git clone https://github.com/samplchallenges/SAMPL-containers.git`
4. Navigate to the "SAMPL-containers/tutorials" directory.
    * command: `cd SAMPL-containers/tutorials`
5. Create a directory called "adv-tut-singularity-base"
    * command: `mkdir adv-tut-singularity-base`
6. Change directories to "adv-tutorial-base"
    * command: `cd adv-tut-singularity-base`

### 1.2: Run a pre-made Docker container and create a conda environment
> In 1.2, we will run the pre-made miniconda container, continuumio/miniconda3, which contains a pre-installed version of miniconda, in interactive mode. This will allow us to interact with the container's command line and directory contents. We will also be able to dynamically create the conda environment we need on the command line inside the container. Even if you have a conda environment installed locally, you will need to complete this step. The container is isolated from your local environment, so it will not have access to your local conda environment.
> 
> Because our container will build off of the miniconda container using it as a base, any environment we create while interatively using the miniconda container should install into our container (which uses miniconda as a base) without additional issues. Building a conda environment outside the miniconda container often results in multiple rounds of trial and error and incompatible packages. We've found the following steps to be the fastest procedure. For more detailed/generalized instructions please see CondaEnvInstructions.pdf.
> 
> When building your own Docker container, this is where you would create your own conda environment with the packages you will need.

1. Start the container. Upon running this command your command line prompt should change. This means you are now inside the container interacting with it's command line and contents. The change in command prompts should look similar to the code block below.
    * command: `singularity shell docker://continuumio/miniconda3`
    ```
    vagrant@ubuntu-bionic:adv-tut-singularity-base$ singularity shell docker://continuumio/miniconda3
    INFO:    Using cached SIF image
    Singularity> 
    ```
2. Create a python 3.7 environmentd
     * command: `conda create -n py37 python=3.7`
3. Activate the environment
     * command: `source activate py37`
4. Install mdtraj and rdkit into the base environment
     * command: `conda install -c conda-forge mdtraj rdkit`
5. Install openbabel 
     * command: `conda install -c openbabel openbabel`
6. Start up the Python interpreter and ensure your version is `3.7.*`. The Python version is 3.7.13 in the code block below.
     ```
     Singularity> python
     Python 3.7.13 (default, Mar 29 2022, 02:18:16) 
     [GCC 7.5.0] :: Anaconda, Inc. on linux
     Type "help", "copyright", "credits" or "license" for more information.
     >>>
     ```
7. In the Python interpreter, import rdkit and mdtraj to ensure there are no errors.
     ```
     >>> import rdkit
     >>> import mdtraj
     ```
8. Quit the Python interpeter
     ```
     >>> quit()
     ```
9. Run openbabel with the help flag to ensure openbabel has installed properly. If properly installed the output should look similar to the code block below
     ```
     Singularity> obabel -H
     Open Babel converts chemical structures from one file format to another

     Usage: 
     obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
     The extension of a file decides the format, unless it is overridden
      by -i or -o options, e.g. -icml, or -o smi
     See below for available format-types, which are the same as the 
     file extensions and are case independent.
     If no input or output file is given stdin or stdout are used instead.
     ```
10. Exit the container shell. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
     * command: `exit`

1.3: Install conda environment (from Section 1.2) into your container
> We will begin creating a `buildfile` to create a container with the virtual environment from the previous section. The only difference in the steps is this time we will update the base environment using `conda update` rather than create a new environment using `conda create`. In the previous step, we could not install into the base environment due to root permissions. 

1. Create and open a file called "buildfile"
2. Copy the following lines from the codeblock into `buildfile`. The following commands contain the instructions to install the conda environment when your container is built
     ```bash
     # Start building off of the docker container continuumio/miniconda3
     Bootstrap: docker
     From: continuumio/miniconda3

     %files

     %post
     # Install python 3.7 and package dependencies into the conda base environment
     conda update conda && \
         conda install python=3.7 && \
         conda install -c conda-forge mdtraj rdkit && \
         conda install -c openbabel openbabel && \
         conda clean --all --yes
     
     %runscript
     ```
     * The `Bootstrap` and `From` lines determine which container to use as a foundation for our container. 
     * The `%files` section states which files to copy into the container ahead of running installation steps. Right now we do not require any files. 
     * The `%post` section is for any installation steps that should occur during the build ahead of runtime. 
     * The `%runscript` section is for any steps that should be executed by the container at runtime.
3. Save `buildfile` and exit your text editor
4. Build your container
     * command: `singularity build --fakeroot adv-tut-base.sif buildfile`
5. Confirm your build succeded by running the container in interactive mode. Upon running this command your command line prompt should change.
     * command: `singularity shell adv-tut-base.sif`
     ```
     vagrant@ubuntu-bionic:adv-tut-singularity-base$ singularity shell adv-tut-base.sif
     INFO:    Using cached SIF image
     Singularity> 
     ```
6. Start up the Python interpreter and ensure your version is `3.7.*`. The Python version is 3.7.13 in the code block below.
     ```
     Singularity> python
     Python 3.7.13 (default, Mar 29 2022, 02:18:16) 
     [GCC 7.5.0] :: Anaconda, Inc. on linux
     Type "help", "copyright", "credits" or "license" for more information.
     >>>
     ```
7. In the Python interpreter, import rdkit and mdtraj to ensure there are no errors.
     ```
     >>> import rdkit
     >>> import mdtraj
     ```
8. Quit the Python interpeter
     ```
     >>> quit()
     ```
9. Run openbabel with the help flag to ensure openbabel has installed properly. If properly installed the output should look similar to the code block below
     ```
     Singularity> obabel -H
     Open Babel converts chemical structures from one file format to another

     Usage: 
     obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
     The extension of a file decides the format, unless it is overridden
      by -i or -o options, e.g. -icml, or -o smi
     See below for available format-types, which are the same as the 
     file extensions and are case independent.
     If no input or output file is given stdin or stdout are used instead.
     ```
10. Exit the container shell. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
     * command: `exit`

  
