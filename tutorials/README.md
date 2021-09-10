# Introduction to SAMPL Containerized Methods

## Purpose:
In [SAMPL4](https://link.springer.com/article/10.1007%2Fs10822-013-9702-2), we learned that human knowledge can be a key factor influencing the success of a computational drug discovery method. To work around this finding, we are creating an automated arm of SAMPL challenges to run methods head-to-head without human intervention. To accomplish this, we will use containerized [Docker](https://www.docker.com/resources/what-container) methods.

The following tutorial is meant to teach the basics of building a simple container to predict docking poses. Here, we use both Python code and command-line programs (specifically, Autodock Vina and MGL Tools).

## Background and Prerequisites

### Important Note on "Docker" versus "Docking":
Please note that "Docker" and "docking" are two separate things.
* **"Docker"** is a program that allows you to containerize methods, essentially allowing you to distribute your method or approach in a reusable way where it can be used reproducibly without human intervention.
* **"Docking"** describes predicting the structure of a complex, in this case a protein-ligand complex.

This terminology is unfortunately not something we can change.

### Expected Background Knowledge
* Basic knowledge of [Python](https://www.python.org/)
* Basic knowledge of [Linux/UNIX command line](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)


### Software Requirements
* Linux or UNIX operating system
* [Docker Desktop](https://www.docker.com/products/docker-desktop)
* [Docker SDK for Python](https://pypi.org/project/docker/)
* [Python 3](https://www.python.org/downloads/)

## Getting started

### Brief Docker Background
* If you have know nothing about Docker, [this video](https://www.youtube.com/watch?v=_dfLOzuIg2o) is a great summary in just 5 minutes.
* Docker containers isolate their internal software from the environment (i.e. someone's operating system or virtual environment) and ensure that the container software works consistently across platforms despite differences in development and staging.
* A [Docker (container) image](https://www.docker.com/resources/what-container) is the blueprint of a Docker container that is not currently running. Docker images contain the instructions to execute your code as a Docker container. 
* The instructions to prepare/build a Docker image are contained in a file called "[Dockerfile](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image)".
* To use a Docker image, we must first build the image with the [`docker build`](https://docs.docker.com/engine/reference/commandline/build/) command and the instructions outlined in Dockerfile. Ensure you are in the proper directory with your Dockerfile and code, then run `docker build -t <name>:<tag/version> .`
   * Examples: `docker build -t adv:0.1 .` or `docker build -t adv:latest .`
* Use the command [`docker images`](https://docs.docker.com/engine/reference/commandline/images/) to list out images you have built or downloaded.
* To delete Docker images, use `docker images` to list your current images and their IMAGE IDs, then run the command [`docker image rm <IMAGE IDs>`](https://docs.docker.com/engine/reference/commandline/image/).

### Pre-Built Autodock Vina Container
A working version of the Autodock Vina container we will build in this tutorial can be found at [Docker Hub under osatom/adv-tutorial](https://hub.docker.com/repository/docker/osatom/adv-tutorial). To play with this container, please use the following steps:
1. Use the [`docker pull`](https://docs.docker.com/engine/reference/commandline/pull/) command to download the "adv-tutorial" Docker container: `docker pull osatom/adv-tutorial:latest`
2. Change directories into the "SAMPL-containers/tutorials" directory: `cd SAMPL-containers/tutorials`
3. Run the command: `python ever_given/run.py osatom/adv-tutorial:latest --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
4. The results will be stored in the directory "tutorials/evergiven_output"


# Tutorial: Build an AutoDock Vina Containerized Method
> This tutorial is separated into five parts: (1) the virtual environment and dependency container build, (2) the docking container build, (3) docking using the Docker container, (4) building your own docking container, and (5) miscellaneous important information. 
> 
> We have separated the base virtual environment build (conda environment and command-line programs) from the implementation build (docking) to illustrate container inheritance and to improve the time required to build the container. In Section 1, we will build an adv-tutorial-base container that includes all the environment and software installations necessary to run AutoDock Vina from a Python script. In Section 2, we will inherit the environment from the Section 1 adv-tutorial-base container to write and build our run-able docking code in the adv-tutorial container. 
>
> This inheritance scheme is also a good practice to improve the time required to build the docking container. When you are writing your docking container, most of the time, you will only make changes to the new docking code you have written and thus only need to re-build the new docking code, rather than re-build the new code AND the environment/software installations. By separating the docking container build from the environment container build, we can build just the docking container each time we change the docking code, improving build time. For example, building only new code (adv-tutorial) typically takes under 5 seconds, but building new code plus environment and software installations (adv-tutorial plus adv-tutorial-base) takes upwards of 115 seconds. This can be analogized to independently compiling parts of a program to improve compile time as a whole later.



## Outline:
* [Section 1: Build the Autodock Vina base container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#section-1-build-the-autodock-vina-base-container)
   * [1.1 Setup](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#11-setup)
   * [1.2 Run a pre-made Docker container and create a conda environment](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#12-run-a-pre-made-docker-container-and-create-a-conda-environment)
   * [1.3 Install conda environment (from Section 1.2) into your container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#13-install-conda-environment-from-section-12-into-your-container)
   * [1.4 Download and prepare the command line programs Autodock Vina and MGL Tools executables for use in the docking container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#14-download-and-prepare-the-command-line-programs-autodock-vina-and-mgl-tools-executables-for-use-in-the-docking-container)
   * [1.5 Install Autodock Vina and MGL Tools into your container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#15-install-autodock-vina-and-mgl-tools-into-your-container)
   * [1.6 Build the base container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#16-build-the-base-container)
* [Section 2: Build the Autodock Vina Docking methods container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#section-2-build-the-container-with-autodock-vina-docking-methods)
   * [2.1 Setup](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#21-setup)
   * [2.2 Add the docking code](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#22-add-the-docking-code)
   * [2.3 Create a setup.py file](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#23-create-a-setuppy-file)
   * [2.4 Write a Dockerfile with instructions to build your container.](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#24-write-a-dockerfile-with-instructions-to-build-your-container)
   * [2.5 Build the docking container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#25-build-the-docking-container)
 * [Section 3: Test/Run your container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#section-3-testrun-your-container)
 * [Section 4: Building Your Own Docking Container](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#section-4-building-your-own-docking-container)
 * [Section 5: Other Important Information](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#section-5-other-important-information)


## Section 1: Build the Autodock Vina base container
> In Section 1, we will build a base container with all necessary packages and programs installed. This base container will serve as a starting point for our docking container. This way, as we write our docking code, it will build quickly since only the docking code will need to be built. 
> 
> For more information, please see [Tutorial: Build an AutoDock Vina Containerized Method](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#tutorial-build-an-autodock-vina-containerized-method)).


### 1.1: Setup
1. Open a terminal
2. Clone the SAMPL-containers repository
    * command: `git clone https://github.com/samplchallenges/SAMPL-containers.git`
3. Navigate to the "SAMPL-containers/tutorials" directory.
4. Create a directory called "adv-tutorial-base" 
   * command: `mkdir adv-tutorial-base`
5. Change directories to "adv-tutorial-base"
   * command: `cd adv-tutorial-base`


### 1.2: Run a pre-made Docker container and create a conda environment

> In 1.2, we will run the pre-made [miniconda](https://docs.conda.io/en/latest/miniconda.html) container, [`continuumio/miniconda3`](https://hub.docker.com/r/continuumio/miniconda3), which contains a pre-installed version of miniconda, in [interactive mode](https://docs.docker.com/engine/reference/run/#foreground). This will allow us to interact with the container's command line and directory contents. We will also be able to dynamically create the conda environment we need on the command line inside the container. **Even if you have a conda environment installed locally, you will need to complete this step.** The container is isolated from your local environment, so it will not have access to your local conda environment. 
> 
> Because our container will build off of the miniconda container using it as a base, any environment we create while interatively using the miniconda container should install into our container (which uses miniconda as a base) without additional issues. Building a conda environment outside the miniconda container often results in multiple rounds of trial and error and incompatible packages. We've found the following steps to be the fastest procedure. For more detailed/generalized instructions please see [CondaEnvInstructions.pdf](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/CondaEnvInstructions.pdf).
>
> When building your own Docker container, this is where you would create your own conda environment with the packages you will need. 
1. Start up Docker Desktop, which will start the Docker daemon. If this is your first time starting Docker Desktop, the application may need to complete some additional installation steps or updates. 
2. Start the container. Upon running this command your command line prompt should change. This means you are now inside the container interacting with it's command line and contents. The change in command prompts should look similar to the code block below. 
   * command: `docker run -it --rm continuumio/miniconda3`
   ```
   megosato@Admins-MacBook-Pro ~ % docker run -it --rm continuumio/miniconda3
   root@7f02be71557e:/# 
   ```
3. Create a conda env called "advenv", if prompted "Proceed ([y]/n)?" type "y"
   * command: `conda create --name advenv python=3.6`
4. Activate advenv: 
   * command: `conda activate advenv`
6. Install [rdkit](https://www.rdkit.org/), if prompted "Proceed ([y]/n)?" type "y"
   * command: `conda install -c conda-forge rdkit`
7. Install [mdtraj](https://mdtraj.org/1.9.4/installation.html), if prompted "Proceed ([y]/n)?" type "y"
   * command: `conda install -c conda-forge mdtraj`
8. Install [openbabel](https://openbabel.org/docs/dev/Installation/install.html), if prompted "Proceed ([y]/n)?" type "y"
   * comand: `conda install -c openbabel openbabel`
9. Start up the Python interpreter and ensure your version is `3.6.*`. The Python version is 3.6.13 in the code block below.
   * command: `python`
   ```
   root@7f02be71557e:/# python
   Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59) 
   [GCC 7.5.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> 
   ```
10. In the Python interpreter, import rdkit, ensure there are no errors.
      ```
      >>> import rdkit
      >>>
      ```
11. In the Python interpreter, import mdtraj, ensure there are no errors.
      ```
      >>> import mdtraj
      >>>
      ```
12. Quit the Python interpeter
      ```
      >>> quit()
      ```
13. Run openbabel with the help flag to ensure openbabel has installed properly. If properly installed the output should look similar to the code block below
      * command: `obabel -H`
      ```
      /opt/app # obabel -H
      Open Babel converts chemical structures from one file format to another

      Usage: 
      obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
      The extension of a file decides the format, unless it is overridden
       by -i or -o options, e.g. -icml, or -o smi
      See below for available format-types, which are the same as the 
      file extensions and are case independent.
      If no input or output file is given stdin or stdout are used instead.
      ```
14. Export the environment. Upon running this command, an organized list of the packages installed into advenv will be printed out to your console. The output will look similar to the condensed output in code block below.
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
15. Copy the output from the export command in step 14 to be pasted into a file in the next section.
17. Exit the container. Upon running this command, you will exit the interactive version of the container and should return to your normal command prompt.
      * command: `exit`


### 1.3: Install conda environment (from [Section 1.2](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#12-starting-a-pre-made-docker-container-and-creating-a-conda-environment)) into your container  

> We will begin creating an environment.yml file that contains all the information about our conda environment packages, and a Dockerfile which contains the instructions required to build the base container. In 1.3, we will only add the necessary commands for installing the conda environment from Section 1.2 to the Dockerfile. We will then test to ensure the container build succeeds at creating the conda environment. 
> 
> For more information on how to write a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image).
1. Create and open a file called "environment.yml" and paste the output you previously copied at Section 1.2, Step 14
2. Change the first line of the file `name: advenv` to `name: base`
3. Delete the last line of the file: `prefix: /opt/conda/envs/advenv`
4. Save the changes to environment.yml and close the file
5. Create and open a file called "Dockerfile"
6. Copy the following lines into Dockerfile. The following commands contain the instructions to install the conda environment when your container is built. 
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

   ENV PATH="/root/.local/bin:$PATH"      
   # set the path
   ```
7. Save the changes to Dockerfile and close the file
8. Build your container into a Docker image to ensure there are no build issues, so far. 
   * command: `docker build -t adv-tutorial-base-test .`
9. If your build from the previous step (step 8) completed without issue, please move on to the next step, otherwise some troubleshooting of the previous steps may be necessary. A successful build looks something like the code block below.
   ```
   (base) megosato@Admins-MacBook-Pro adv-tutorial-base % docker build -t adv-tutorial-base-test .
   [+] Building 2.0s (13/13) FINISHED                                                                                                         
    => [internal] load build definition from Dockerfile                                                                                  0.0s
    => => transferring dockerfile: 402B                                                                                                  0.0s
    => [internal] load .dockerignore                                                                                                     0.0s
    => => transferring context: 2B                                                                                                       0.0s
    => [internal] load metadata for docker.io/continuumio/miniconda3:4.9.2-alpine                                                        1.1s
    => [auth] continuumio/miniconda3:pull token for registry-1.docker.io                                                                 0.0s
    => [internal] load build context                                                                                                     0.8s
    => => transferring context: 74.91MB                                                                                                  0.8s
    => [1/7] FROM docker.io/continuumio/miniconda3:4.9.2-alpine@sha256:82bd96b0e95188e152d137f6c9834ea731bfc78e5c5f27b3c90f2be31e9e61d8  0.0s
    => CACHED [2/7] WORKDIR /opt/app/                                                                                                    0.0s
    => CACHED [3/7] COPY . ./                                                                                                            0.0s
    => CACHED [4/7] RUN conda env update -f environment.yml &&     conda clean --all --yes                                               0.0s
    => CACHED [5/7] RUN pip install .                                                                                                    0.0s
    => CACHED [6/7] RUN /opt/app/dependencies/mgl/install.sh                                                                             0.0s
    => CACHED [7/7] RUN /opt/app/dependencies/adv/bin/vina --help                                                                        0.0s
    => exporting to image                                                                                                                0.0s
    => => exporting layers                                                                                                               0.0s
    => => writing image sha256:0bbc54e5f8a5fbc5274a29fb7e20326309eb01aa2ece18bea1d7406a6fcbecea                                          0.0s
    => => naming to docker.io/library/adv-tutorial-base-test                                                                             0.0s

   Use 'docker scan' to run Snyk tests against images to find vulnerabilities and learn how to fix them
   ```
10. Delete the Docker image previously built in part Step 8. Step 8 was just to ensure nothing had gone wrong in the build so far.
    * List the Docker images: `docker images`
      ```
      REPOSITORY                      TAG       IMAGE ID       CREATED          SIZE
      adv-tutorial-base-test          latest    5cf75f044da5   15 minutes ago   1.04GB
      ```
    * Determine the "IMAGE ID" of "adv-tutorial-base-test", in this case the ID was:  "5cf75f044da5"
    * Delete the Docker image: `docker image rm [IMAGE ID]`

### 1.4: Download and prepare the command line programs Autodock Vina and MGL Tools executables for use in the docking container

> In 1.4, we will incorporate the command line tools [Autodock Vina](http://vina.scripps.edu/) and [MGL Tools](http://mgltools.scripps.edu/downloads) into our base container. Please do not change which installers you download based on your native operating system (OS) because the OS used inside Docker container, Linux x86, may differ from your native OS. For example, I am currently on a Mac, but the OS inside the Docker container is Linux x86, so any installers that work for my native Mac OS would not work inside my Docker container. 
> 
> When building your own container, this is where you would add in any command line program files.

* NOTE: Please note the commands with tarball files may be downloaded as `.tgz`, `.tar`, or `.tar.gz`, please alter the commands according to the name of the file downloaded.
1. Create a directory called "dependencies". Upon creating the dependencies directory, your directory structure should look like the code block below.
   * command: `mkdir dependencies`
   ```
   adv-tutorial-base
   ├── Dockerfile
   ├── dependencies
   └── environment.yml
   ```
2. Download Autodock Tools linux x86 "autodock_vina_1_1_2_linux_x86.tgz" from http://vina.scripps.edu/download.html
3. Move "autodock_vina_1_1_2_linux_x86.tgz" into "adv-tutorial-base"
   * command: `mv {path_to_download}/autodock_vina_1_1_2_linux_x86.tgz .`
5. Untar "autodock_vina_1_1_2_linux_x86.tgz"
   * command: `tar -xvf autodock_vina_1_1_2_linux_x86.tgz`
6. Delete the .tgz file:
   * command: `rm autodock_vina_1_1_2_linux_x86.tgz`
7. Rename "autodock_vina_1_1_2_linux_x86" to "adv"
   * command: `mv autodock_vina_1_1_2_linux_x86 adv`
8. Move "adv" directory to inside "dependencies" directory
   * command: `mv adv dependencies`
9. Download MGL Tools linux x86 `mgltools_x86_64Linux2_1.5.6.tar.gz` from http://mgltools.scripps.edu/downloads
10. Move "mgltools_x86_64Linux2_1.5.6.tar.gz" into "adv-tutorial-base"
      * command: `mv {path_to_download}/mgltools_x86_64Linux2_1.5.6.tar.gz .`
11. Untar "mgltools_x86_64Linux2_1.5.6.tar.gz"
      * command `tar -xvf mgltools_x86_64Linux2_1.5.6.tar.gz`
12. Delete the .tgz file:
      * command: `rm mgltools_x86_64Linux2_1.5.6.tar.gz`
13. Rename "mgltools_x86_64Linux2_1.5.6" to "mgl"
      * command: `mv mgltools_x86_64Linux2_1.5.6 mgl`
14. Move "mgl" directory to inside the "dependencies" directory
      * command: `mv mgl dependencies`
15. Open "dependencies/mgl/install.sh"
16. Change line 6 from `TarDir='pwd'` to `TarDir="/opt/app/dependencies/mgl/"`
17. Change line 7 from `export MGL_ROOT=""` to `export MGL_ROOT="/opt/app/dependencies/mgl/"`
18. Save and close "dependencies/mgl/install.sh"


### 1.5: Install Autodock Vina and MGL Tools into your container
> In 1.5, we will add in the installation commands for Autodock Vina and MGL Tools to the Dockerfile. 
> 
> When building your own container, this is where you would add in any command line program installation steps.
1. Re-open the Dockerfile from [Section 1.3](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#13-install-conda-environment-from-section-12-into-your-container) and append the following lines to the end of the file
   ```
   RUN /opt/app/dependencies/mgl/install.sh
   # command to install MGL Tools

   ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"
   # add MGL Tools and Autodock Vina to your path

   RUN /opt/app/dependencies/adv/bin/vina --help
   # Run autodock vina to ensure its properly installed
   ```
2. Save the changes to the `Dockerfile` and close the file


### 1.6: Build the base container

> In 1.6, we will build a Docker image using "docker build", which we will inherit from in Section 2.

1. Build the base container
   * command: `docker build -t adv-tutorial-base .`
2. Ensure your build was successful. A successful build will look something like the code block below.
   ```
   (base) megosato@Admins-MacBook-Pro adv-tutorial-base % docker build -t adv-tutorial-base .
   [+] Building 117.9s (11/11) FINISHED                                                                            
    => [internal] load build definition from Dockerfile                                                       0.0s
    => => transferring dockerfile: 37B                                                                        0.0s
    => [internal] load .dockerignore                                                                          0.0s
    => => transferring context: 2B                                                                            0.0s
    => [internal] load metadata for docker.io/continuumio/miniconda3:4.9.2-alpine                             0.5s
    => [1/6] FROM docker.io/continuumio/miniconda3:4.9.2-alpine@sha256:82bd96b0e95188e152d137f6c9834ea731bfc  0.0s
    => [internal] load build context                                                                          0.0s
    => => transferring context: 2.62kB                                                                        0.0s
    => CACHED [2/6] WORKDIR /opt/app/                                                                         0.0s
    => CACHED [3/6] COPY . ./                                                                                 0.0s
    => [4/6] RUN conda env update -f environment.yml &&     conda clean --all --yes                         103.0s
    => [5/6] RUN /opt/app/dependencies/mgl/install.sh                                                         3.2s 
    => [6/6] RUN /opt/app/dependencies/adv/bin/vina --help                                                    0.2s 
    => exporting to image                                                                                    10.8s 
    => => exporting layers                                                                                   10.8s 
    => => writing image sha256:594ea10a0d26b3467246416fb1ce93e9c4b4cbb765038038e2020c828fb6c210               0.0s 
    => => naming to docker.io/library/adv-tutorial-base                                                       0.0s 

   Use 'docker scan' to run Snyk tests against images to find vulnerabilities and learn how to fix them 
   ```


## Section 2: Build the container with Autodock Vina Docking methods

> Previously, we built the adv-tutorial-base container which contains all the environment and software installations necessary to run AutoDock Vina from a Python script. Please see [Tutorial: Build an AutoDock Vina Containerized Method](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#tutorial-build-an-autodock-vina-containerized-method) for a more detailed explanation as to why we are using 2 separate builds. In Section 2, we will inherit from the [Section 1](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#section-1-build-the-autodock-vina-base-container) base container to write and build our run-able docking code. 

### 2.1: Setup

1. Navigate up a directory to "tutorials"
   * command: `cd ..`
3. Create a directory called "adv-tutorial" in the tutorials directory
   * command: `mkdir adv-tutorial`
4. Change directories to "adv-tutorial"
   * command: `cd adv-tutorial`

### 2.2: Add the docking code

> In 2.2, we will incorporate the docking code into our container directory. For the sake of simplicity, we will be using pre-written docking code. Please see [tutorials/BuildYourOwnDockingContainer.md](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/BuildYourOwnDockingContainer.md) for more information on the inputs and kwargs required of each main function.
> 
> When building your own container, this is where you would add in your methods.
1. Copy the AutoDock class file from "SAMPL-containers/docking/examples/adv-tutorial/autodock.py" to "adv-tutorial"
   * command: `cp ../../docking/examples/adv-tutorial/autodock.py .`
2. Copy the AutoDock main file from "SAMPL-containers/docking/examples/adv-tutorial/main.py" to "adv-tutorial"
   * command: `cp ../../docking/examples/adv-tutorial/main.py .`

### 2.3: Create a setup.py file

> In 2.3, we will customize a setup.py file to to match the Python modules we have written.
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
2. Modify the "py_modules" list, by adding the two modules ([autodock.py](https://github.com/samplchallenges/SAMPL-containers/blob/main/docking/examples/adv-tutorial/autodock.py) and [main.py](https://github.com/samplchallenges/SAMPL-containers/blob/main/docking/examples/adv-tutorial/main.py)) with our docking code from the previous subsection ([2.2](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#22-add-the-docking-code)): "autodock" and "main". This is where you 
    ```
    py_modules=[
        'autodock',
        'main',
    ]
    ```
3. Modify the "install_requires" list, adding "Click" a pip installable package we did not add to the previous base build
    ```
    install_requires=[
        'Click',
    ]
    ```
4. Modify the "entry_points" section, by adding an [entry point](https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html#console-scripts) in the format `{command-to-call-in-Dockerfile}={py_module_with_main}:{function_to_run}`. The Python module "main.py" contains the main function "main_function".
    ```
    entry_points='''
        [console_scripts]
        run-autodock=main:main_function
    '''
    ```
5. Save and close setup.py

### 2.4: Write a Dockerfile with instructions to build your container

> In 2.4, we will create a Dockerfile which contains the instructions required to build the docking container, as well as the [ENTRYPOINT](https://docs.docker.com/engine/reference/builder/#entrypoint) (see [2.3](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#23-create-a-setuppy-file)) which tells the container which file and function to run. For more information on how to write a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image). 
1. Create and open a file called "Dockerfile"
2. Copy and paste the following into Dockerfile
    ```
    FROM

    WORKDIR /opt/app/
    # set the work directory

    COPY

    RUN pip install .
    # install setup.py

    ENV PATH="/root/.local/bin:$PATH"

    RUN ls -l /opt/app

    ENTRYPOINT []
    ```

2. Next to "FROM" add the name of the base build from Section 1. This is where we specify that our container will inherit from the base we built in [Section 1](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials#section-1-build-the-autodock-vina-base-container)
   * `FROM adv-tutorial-base`
3. Next to "COPY" add the names of all files necessary to run our docking program, including "setup.py", "autodock.py", and "main.py" into the container directory  "./" or "/opt/app". This is where you would specify the files you have wirtten for your docking program in place of "autodock.py" and "main.py"
   * `COPY setup.py autodock.py main.py ./`
4. Next to "ENTRYPOINT", add the "entry_point" you declared in step 4 of the previous subsection ([2.3](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#23-create-a-setuppy-file)) inside the brackets, in quotations.
   * `ENTRYPOINT ["run-autodock"]`
5. Save and close Dockerfile


### 2.5: Build the docking container

> In 2.5, we will build a Docker image that will execute our docking program when run. For more information about the `docker build` command, please see the [official Docker documentation](https://docs.docker.com/engine/reference/commandline/build/).
1. Build the container. A successful build should look similar to the code block below.
   * command: `docker build -t adv-tutorial .`
   ```
   (base) megosato@Admins-MacBook-Pro adv-tutorial % docker build -t adv-tutorial .
   [+] Building 1.6s (10/10) FINISHED                                                                                                         
    => [internal] load build definition from Dockerfile                                                                                  0.0s
    => => transferring dockerfile: 227B                                                                                                  0.0s
    => [internal] load .dockerignore                                                                                                     0.0s
    => => transferring context: 2B                                                                                                       0.0s
    => [internal] load metadata for docker.io/osatom/adv-rdkit-base:latest                                                               1.5s
    => [auth] osatom/adv-rdkit-base:pull token for registry-1.docker.io                                                                  0.0s
    => [1/4] FROM docker.io/osatom/adv-rdkit-base:latest@sha256:2845b4f7491dc5df3777c70a4ac2643a093b49b424c7c9f3748b1b85e17bfa72         0.0s
    => => resolve docker.io/osatom/adv-rdkit-base:latest@sha256:2845b4f7491dc5df3777c70a4ac2643a093b49b424c7c9f3748b1b85e17bfa72         0.0s
    => [internal] load build context                                                                                                     0.0s
    => => transferring context: 11.98kB                                                                                                  0.0s
    => CACHED [2/4] WORKDIR /opt/app/                                                                                                    0.0s
    => CACHED [3/4] COPY setup.py autodock.py main.py ./                                                                                 0.0s
    => CACHED [4/4] RUN pip install .                                                                                                    0.0s
    => exporting to image                                                                                                                0.0s
    => => exporting layers                                                                                                               0.0s
    => => writing image sha256:6d41aab6331068fd765f47ef1d48467567f62595e52c629a4caeadb57aa7740c                                          0.0s
    => => naming to docker.io/library/adv-tutorial                                                                                       0.0s

   Use 'docker scan' to run Snyk tests against images to find vulnerabilities and learn how to fix them
   ```

## Section 3: Test/Run your container
In this section, we will use the wrapper `ever_given` to run the docking container. `ever_given` mimics the infrastructure we will use to run your container on the SAMPL-league website, making it a great way to test that you container will run properly ahead of uploading to the [SAMPL challenges website](https://app.samplchallenges.org/). `ever_given` also abstracts away [volume mounting](https://docs.docker.com/engine/reference/commandline/run/#mount-volume--v---read-only) to link your local directory with a directory inside the container, making it easier to quickly test your container. For more information on how to use `ever_given` please see [ever_givenUsage.md](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/ever_givenUsage.md).
1. Change directories into "tutorials" one directory above:
   * command: `cd ..`
2. Run the container. A successful run should look like the code block below.
   * command: `python ever_given/run.py adv-tutorial --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
   ```
   (base) megosato@Admins-MacBook-Pro tutorials % python ever_given/run.py adv-tutorial --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 
   --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor
   command? 
   file kwargs {'receptor': 'data/receptor.pdb', 'hint': 'data/hint.pdb'}
   kwargs {'hint_radius': '6', 'hint_molinfo': 'E51', 'smiles': 'CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O'}
   Putting output into evergiven_output
   evergiven_output docked_ligand /mnt/outputs/best_dock.pdb
   evergiven_output receptor /mnt/outputs/rec-dock.pdb
   Results: {'docked_ligand': 'evergiven_output/best_dock.pdb', 'receptor': 'evergiven_output/rec-dock.pdb'}
   ```
3. Examine your container outputs in the directory "tutorials/evergiven_output"
* Note: If you receive an error similar to the one below, please ensure [Docker SDK](https://pypi.org/project/docker/) is installed. 
    ```
    Traceback (most recent call last):
    File "/Users/megosato/SAMPL-containers/tutorials/ever_given/run.py", line 4, in <module>
      from ever_given import wrapper
    File "/Users/megosato/SAMPL-containers/tutorials/ever_given/ever_given/wrapper.py", line 8, in <module>
      import docker
    ModuleNotFoundError: No module named 'docker'
    ```


## Section 4: Building Your Own Docking Container
* For more detailed information about docking container requirements and how to modify this tutorial to suit your needs, please see [BuildYourOwnDockingContainer.md](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/BuildYourOwnDockingContainer.md)
* For an example template of a docking container directory, please see [SAMPL-containers/tutorial/templates/docking](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials/templates/docking)
* For more detailed information about LogD container requirements and how to modify this tutorial to suit your needs [BuildYourOwnLogDContainer.md](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/BuildYourOwnLogDContainer.md)
* For an example template of a LogD container directory, please see [SAMPL-containers/tutorial/templates/logd](https://github.com/samplchallenges/SAMPL-containers/tree/main/tutorials/templates/logd)
* For more information on how to build your own conda environment inside a container, please see [CondaEnvInstructions.pdf](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/CondaEnvInstructions.pdf)


## Section 5: Other Important Information
* For more detailed tutorials on how to use Docker please see the following resources:
  * [Official Docker Documentation](https://docs.docker.com/get-started/)
  * [Brief Docker Tutorial (12m)](https://www.youtube.com/watch?v=YFl2mCHdv24)
  * [Docker Beginners Course (2hrs)](https://www.youtube.com/watch?v=fqMOX6JJhGo)
  * [Brief video on how Docker containers work (15m)](https://www.youtube.com/watch?v=rOTqprHv1YE)
  * YouTube is a great resource for learning Docker, feel free to search for other tutorials that suit your specific needs as well
* To upload your own container to docker.io create a [Docker Hub account](https://hub.docker.com/signup)
  * Please note that free Docker Hub accounts only receive 1 private repository.
  * If you have sensitive information inside your containers, please ensure you are using a private repository.
* Eventually, upon building enough Docker images, you may begin to run out of memory. Please remember to regularly delete any Docker images you no longer need. (See [Brief Docker Usage Tips](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/README.md#brief-docker-background): last bullet point)
* If you get an error similar to `Error response from daemon: Bad response from Docker engine` when using a Docker command, make sure your Docker daemon is running by starting Docker Desktop. 
