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
        * `'singularity_local'`: sif file image (ending in `.sif` built with singularity engine
    * `{file-kwargs}` are keyword option-value pairs where the value is a file path the container will need access to. Any keyword that expects a file path should get `--file-` prepended to it. 
    * `{kwargs}` are keyword option-value pairs where the value is anything other than a file path.
    * `--output-keys {keys}` are any output keys as a comma separated list from required the print output statements that contain file paths. `ever_given` automatically handles the file path mounting and unmounting, and will update the file path for us
