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
* A Docker image contains the instructions to execute your code inside a docker container. Docker images will store the information necessary to run your containerized methods
* To use an image, we must first build it. Ensure you are in the directory with your Dockerfile and container code, then run `docker build -t <name>:<tag/version> .`
   * Examples: `docker build -t adv:0.1 .` or `docker build -t adv:latest .`
* Use the command docker images to list out your built images
* To delete Docker images, use docker images to list your current images and their IMAGE IDs, then run the command `docker image rm <IMAGE IDs>`

### Pre-Built Autodock Vina Container
A working version of the Autodock Vina container we will build in this tutorial can be found at [Docker Hub under osatom/adv-tutorial](https://hub.docker.com/repository/docker/osatom/adv-tutorial). To play with this container, please use the following steps:
1. Pull the "adv-tutorial" docker container: `docker pull osatom/adv-tutorial:latest`
2. Change directories into the "examples" directory: `cd examples`
3. Run the command: `python ever_given/run.py osatom/adv-tutorial:latest --file-receptor data/receptor.pdb --file-hint data/hint.pdb --hint_radius 6 --hint_molinfo "E51" --smiles "CCCCNc1cc(cc(n1)OC)C(=O)N[C@@H](Cc2ccccc2)[C@H](C[C@@H](C)C(=O)NCCCC)O" --output-keys docked_ligand,receptor`
4. The results will be stored in the directory "examples/evergiven_output"


# Tutorial: Build an AutoDock Vina Containerized Method
> In this tutorial, we have separated the virtual environment build (i.e. conda environment) from the implementation build (docking code) to illustrate container inheritance and improve the time required to build the container. In Section 1, we will build an adv-tutorial-base container which contains all the environment and software installations necessary to run AutoDock Vina from a python script.  In Section 2, we will inherit from the Section 1 adv-tutorial-base container to write and build our run-able docking code in the adv-tutorial container. This inheritance scheme is also a good practice to improve time required to build the container. When you are writing, modifying, and testing your own docking container, each time you build your container to test it, you will only need to build the new code you have written, rather than the new code AND the environment/software installations. For example, for the AutoDock Vina Tutorial container building just the new code (adv-tutorial) typically takes under 5 seconds, but building the new code plus the environment and software installations (adv-tutorial plus adv-tutorial-base) takes upwards of 115 seconds. This can be analogized to independently compiling parts of a C program to improve build time later.


## Outline:
* [Section 1: Build the Autodock Vina base container](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#section-1-build-the-autodock-vina-base-container)
   * [1.1 Setup](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#11-setup)
   * [1.2 Starting a docker file and creating a conda environment](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#12-create-a-conda-environment)
   * [1.3 Install conda environment (from Section 1.2) into your container](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/README.md#13-install-conda-environment-from-section-12-into-your-container)
   * [1.4 Install the command line programs: Autodock Vina and MGL Tools executables](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/README.md#14-install-the-command-line-programs-autodock-vina-and-mgl-tools-executables)
   * [1.5 Update the Dockerfile to include Autodock Vina and MGL Tools installation steps](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#15-update-the-dockerfile-to-include-autodock-vina-and-mgl-tools-installations)
   * [1.6 Build the base container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#16-build-the-base-container)
* [Section 2: Build the Autodock Vina Docing methods container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#section-2-build-the-container-with-autodock-vina-docking-methods)
   * [2.1 Setup](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#21-setup)
   * [2.2 Add the docking code](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#22-add-the-docking-code)
   * [2.3 Create a setup.py file](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#23-create-a-setuppy-file)
   * [2.4 Write a Dockerfile with instructions to build your container.](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/README.md#24-write-a-dockerfile-with-instructions-to-build-your-container)
   * [2.5 Build the docking container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#24-create-a-dockerfile)
 * [Section 3: Test/Run your container](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#section-3-testrun-your-container)


## Section 1: Build the Autodock Vina base container
> In this section, we will build a base container that has all necessary packages and programs installed. This way, as we write our docking code it will build quickly since only the docking code will need to be built.


### 1.1: Setup
1. Open a terminal and navigate to the "SAMPL-containers/tutorial" directory.
2. Create a directory called "adv-tutorial-base" 
   * command: `mkdir adv-tutorial-base`
2. Change directories to "adv-tutorial-base"
   * command: `cd adv-tutorial-base`


### 1.2: Run a pre-made docker container and create a conda environment

> In 1.2, we will run the pre-made [miniconda](https://docs.conda.io/en/latest/miniconda.html) container, "continuumio/miniconda3", which contains a pre-installed version of miniconda, in [interactive mode](https://docs.docker.com/engine/reference/run/#foreground). This will allow us to interact with the container's command line and directory contents. We will also be able to dynamically create the conda environment we need on the command line inside the container. Even if you have a conda environment installed locally, you will need to complete this step. The container is an isolated virtual machine, so it will not have access to your local conda environment. 
> 
> Because our container will build off of the miniconda container using it as a base, any environment we create while interatively using the miniconda container should install into our container without additional issues. Building a conda environment outside the miniconda container often results in multiple rounds of trial and error and incompatible packages. We've found the following steps to be the fastest procedure. For more detailed/generalized instructions please see [CondaEnvInstructions.pdf](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/CondaEnvInstructions.pdf)
>
> When building your own docker container, this is where you would create your own conda environment with the packages you will need. 
1. Start up Docker Desktop, which will start the Docker daemon. If this is your first time starting Docker Desktop, the application may need to complete some additional installation steps or updates. 
2. Start the container. Upon running this command your command line prompt should change. This means you are now inside the container interacting with it's command line and contents. 
   * command: `docker run -it --rm continuumio/miniconda3`
   ```
   megosato@Admins-MacBook-Pro ~ % docker run -it --rm continuumio/miniconda3
   root@7f02be71557e:/# 
   ```
3. Create a conda env called "advenv", if prompted "Proceed ([y]/n)?" please type "y"
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
14. Export the environment. Upon running this command, an organized list of the environment packages will be printed out to your console (something like the code block below).
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


### 1.3: Install conda environment (from [Section 1.2](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#12-starting-a-pre-made-docker-container-and-creating-a-conda-environment)) into your container  

> We will begin creating an environment.yml file which contains all the information about our conda environment packages, as well as a Dockerfile which contains the instructions required to build the base container. In 1.3, we will only add the necessary commands for installing the conda environment from Section 1.2 to the Dockerfile. We will then test to ensure the build with the conda environment installed succeeds. For more information on how to write a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image).
1. Create and open a file called "environment.yml" and paste the output you previously copied at Section 1.2 Step 15
2. Change the first line of the file `name: advenv` to `name: base`
3. Delete the last line of the file: `prefix: /opt/conda/envs/advenv`
4. Save the changes to environment.yml and exit
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
3. Save the changes to Dockerfile and exit
4. Build your container into a docker image to ensure there are no build issues, so far. 
   * command: `docker build -t adv-tutorial-base-test .`
5. If your build from the previous step (step 4) completed without issue, please move on to the next step, otherwise some troubleshooting may be necessary. A successful build looks something like the code block below.
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
6. Delete the docker image previously built in part Step 4. Step 4 was just to ensure nothing has gone wrong in the build so far.
   * List the docker images: `docker images`
      ```
      REPOSITORY                      TAG       IMAGE ID       CREATED          SIZE
      adv-tutorial-base-test          latest    5cf75f044da5   15 minutes ago   1.04GB
      ```
   * Determine the "IMAGE ID" of "adv-tutorial-base-test", in my case it was:  "5cf75f044da5"
   * Delete the docker image: `docker image rm [IMAGE ID]`

### 1.4: Install the command line programs: Autodock Vina and MGL Tools executables

> In 1.4, we will incorporate the command line tools [Autodock Vina](http://vina.scripps.edu/) and [MGL Tools](http://mgltools.scripps.edu/downloads) into our base container. Please do not change which installers you download based on your native operating system (OS) because the OS used inside docker container, Linux x86, may be different than your native OS. For example, I am currently on a Mac, but the OS inside the docker container is Linux x86, so any installers that work for my native Mac OS would not work inside my docker container. 
> 
> When building your own container, this is where you would add in any command line program files.
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
18. Close and save "dependencies/mgl/install.sh"


### 1.5: Update the Dockerfile to include Autodock Vina and MGL Tools installation steps
> In 1.5, we will add in the installation commands for Autodock Vina and MGL Tools to the Dockerfile. 
> 
> When building your own container, this is where you would add in any command line program installation steps.
1. Re-open the Dockerfile from [Section 1.3](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#13-create-a-dockerfile) and append the following lines to the end of the file
   ```
   RUN /opt/app/dependencies/mgl/install.sh
   # command to install MGL Tools

   ENV PATH="/opt/app/dependencies/:/opt/app/dependencies/mgl/:$PATH"
   # add MGL Tools and Autodock Vina to your path

   RUN /opt/app/dependencies/adv/bin/vina --help
   # Run autodock vina to ensure its properly installed
   ```
2. Save the changes to the `Dockerfile` and exit


### 1.6: Build the base container

> In 1.6, we will build a docker image using "docker build", which we will inherit from in the Section 2.

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

> Previously, we built the adv-tutorial-base container which contains all the environment and software installations necessary to run AutoDock Vina from a python script. Please see [Tutorial: Build an AutoDock Vina Containerized Method](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#tutorial-build-an-autodock-vina-containerized-method) for a more detailed explanation as to why we use 2 separate builds. In Section 2, we will inherit from the [Section 1](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#section-1-build-the-autodock-vina-base-container) base container to write and build our run-able docking code. 

### 2.1: Setup

1. Navigate up a directory to "examples"
   * command: `cd ..`
3. Create a directory called "adv-tutorial" in the examples directory
   * command: `mkdir adv-tutorial`
4. Change directories to "adv-tutorial"
   * command: `cd adv-tutorial`

### 2.2: Add the docking code

> In 2.2, we will incorporate the docking code into our container directory. For the sake of simplicity, we will be using pre-written docking code. Please see [examples/ContainerRequirements.md](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/ContainerRequirements.md) for more information on the inputs and kwargs required of each main function.
> 
> When building your own container, this is where you would add in your methods.
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
2. Modify the "py_modules" list, by adding the two modules (autodock.py and main.py) with our docking code from the previous subsection ([2.2](https://github.com/samplchallenges/SAMPL-league/tree/containers/examples#22-add-the-docking-code)): "autodock" and "main". This is where you 
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
4. Modify the "entry_points" section, by adding an [entry point](https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html#console-scripts) in the format `{command-to-call-in-Dockerfile}={py_module_with_main}:{function_to_run}`. The python module "main.py" contains the main function "main_function".
    ```
    entry_points='''
        [console_scripts]
        run-autodock=main:main_function
    '''
    ```
5. Save and close setup.py

### 2.4: Write a Dockerfile with instructions to build your container

> In 2.4, we will create a Dockerfile which contains the instructions required to build the docking container, as well as the [ENTRYPOINT](https://docs.docker.com/engine/reference/builder/#entrypoint) (see [2.3](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#23-create-a-setuppy-file)) which tells the container which file and function to run. For more information on how to visit a Dockerfile, please see the [official Docker documentation](https://docs.docker.com/get-started/02_our_app/#build-the-apps-container-image). 
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

2. Next to "FROM" add the name of the base build from Section 1. This is where we specify that our container will inherit from the base we built in [Section 1](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials#section-1-build-the-autodock-vina-base-container)
   * `FROM adv-tutorial-base`
3. Next to "COPY" add the names of all files necessary to run our docking program, including "setup.py", "autodock.py", and "main.py" into the container directory  "./" or "/opt/app". This is where you would specify the files you have wirtten for your docking program in place of "autodock.py" and "main.py"
   * `COPY setup.py autodock.py main.py ./`
4. Next to "ENTRYPOINT", add the "entry_point" you declared in step 4 of the previous subsection ([2.3](https://github.com/samplchallenges/SAMPL-league/blob/containers/examples/README.md#23-create-a-setuppy-file)) inside the brackets, in quotations.
   * `ENTRYPOINT ["run-autodock"]`
5. Save and close Dockerfile


### 2.5: Build the docking container

> In 2.5, we will build a docker image that will execute our docking program when run. For more information about the `docker build`, please see the [official Docker documentation](https://docs.docker.com/engine/reference/commandline/build/).
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
In this section, we will use the wrapper `ever_given` to run the docking container. `ever_given` mimics the infrastructure we will use to run your container on the SAMPL-league website, making it a great way to test that you container will run properly ahead of uploading to the [SAMPL challenges website](https://app.samplchallenges.org/). `ever_given` also abstracts away [volume mounting](https://docs.docker.com/engine/reference/commandline/run/#mount-volume--v---read-only) to link your local directory with a directory inside the container, making it easier to quickly test your container. 
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


# Building Your Own Docking Container
* For more detailed information about container requirements, please see [ContainerRequirements.md](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/ContainerRequirements.md)
* For an example template of a container directory, please see [SAMPL-containers/tutorial/template/](https://github.com/samplchallenges/SAMPL-containers/tree/tutorial/tutorials/template)
* For more information on how to build your own conda environment inside a container, please see [CondaEnvInstructions.pdf](https://github.com/samplchallenges/SAMPL-containers/blob/tutorial/tutorials/CondaEnvInstructions.pdf)


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
