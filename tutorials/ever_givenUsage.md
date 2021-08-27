# `ever_given` Usage Notes

`ever_given` is a simple wrapper to run a docker container. It abstracts away dealing with [volume mounting](https://docs.docker.com/storage/volumes/) allowing files used or made by the container persist beyond the container's runtime when necessary.

## Running a container with ever_given
* To run a container with ever_given use: `python run.py {container_name}:{version/tag} {file-kwargs} {kwargs} --output-keys {keys}`
    * `{file-kwargs}` are keyword option, value pairs where the value is a file path the container will need access to. Any keyword that expects a file path should get `--file-` prepended to it. 
        * Example: In the Docking container example, there are 2 keyword that are files `--hint` and `--receptor`. When running with `ever_given`, these would become `--file-hint` and `--file-receptor`
    * `{kwargs}` are keyword option, value pairs where the value is anything other than a file path.
         * Example: In the Docking container example, `--hint_radius`, `--hint_molinfo` and `--smiles` are all passed as is with no changes to the keyword.
    * `--output-keys {keys} are any output keys as a comma separated list from required the print output statements that contain file paths. `ever_given` automatically handles the file path mounting and unmounting, and will update the file path for us
         * Example: In the Docking container example, both output keys, `docked_ligand` and `receptor`, are paired with file paths. When running with `ever_given we would add `--output_keys docked_ligand,receptor` to our `python run.py` command
