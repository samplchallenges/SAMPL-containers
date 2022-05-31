# `ever_given` Usage Notes

`ever_given` is a simple wrapper to run docker or singularity containers. It abstracts away [volume mounting](https://docs.docker.com/storage/volumes/) allowing files used or made by the container to persist beyond the container's runtime when necessary with little intervention on your part.

## Installing `ever_given`
* `ever_given` is pip installable
   * `pip install ever-given`

## Running a container with `ever_given`
* To run a container with `ever_given` use:  
    ```
    evergiven {container_name}:{version/tag} --container-engine {engine} --container-type {type} {file-kwargs} {kwargs} --output-keys {keys}
    ```
    * `--container-engine` is the container engine to use.
      * Currently supported:
        * `'docker'`: docker engine can only run docker container images
        * `'singularity'`: singularity engine can run singularity or docker container images
    * `--container-type` is the container type of the image you are running.
      * Currently supported:
        * `'docker'`: container image built with docker engine
        * `'singularity_local'`: sif file image (ending in `.sif`) built with singularity engine
    * `{file-kwargs}` are keyword option-value pairs where the value is a file path the container will need access to. Any keyword that expects a file path should get `--file-` prepended to it. The `--file-` prepend is very important as it allows us
      * Example:
        * `--file-receptor '/path/to/receptor.pdb'`
    * `{kwargs}` are keyword option-value pairs where the value is anything other than a file path.
      * Example:
        * `--smiles 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'`
    * `--output-keys {comma_separated_keys}` are any output keys as a comma separated list from the required print output statements that represent file paths. `ever_given` automatically handles the file path mounting and unmounting
      * Example:
        * `--output-keys file_output_key1,file_output_key2,file_output_key3`
