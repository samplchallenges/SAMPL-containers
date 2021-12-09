# Conda Environment Instructions



## Creating a Conda Environment Inside a Container
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
9. In the directory housing your container Dockerfile, code, etc. create a file called `environment.yml` 
10. Open `environment.yml` for editing and paste in the output from step 6.
11. Change the first line of the file from `name: env_name` to `name: base`
12. Delete the last line of the file: `prefix: /opt/conda/envs/advenv`
  * The file should now look like the figure below.
  ![yaml](https://github.com/samplchallenges/SAMPL-containers/blob/main/tutorials/images/conda_env_yml.png)
14. Save and close `environment.yml`




