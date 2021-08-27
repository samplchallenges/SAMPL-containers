# `ever_given` Usage Notes

`ever_given` is a simple wrapper to run a docker container. It abstracts away dealing with [volume mounting](https://docs.docker.com/storage/volumes/) allowing files used or made by the container persist beyond the container's runtime when necessary.

## Running a container with ever_given
* To run a container with ever_given use: `python run.py {container_name}:{version/tag} {file-kwargs} {kwargs} --output-keys {keys}
    * `{file-kwargs}` are key, value pairs where the value is a file path the container will need access to. 
        * EX: In the Docking container example, there are 2 keyword arguments that are files `--hint` and `--receptor`
