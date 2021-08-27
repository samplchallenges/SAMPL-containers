# `ever_given` Usage Notes

`ever_given` is a simple wrapper to run a docker container. It abstracts away dealing with [volume mounting](https://docs.docker.com/storage/volumes/) allowing files used or made by the container persist beyond the container's runtime when necessary.

## Running a container with `ever_given`
* To run a container with `ever_given` use: `python run.py {container_name}:{version/tag} {file-kwargs} {kwargs} --output-keys {keys}`
    * `{file-kwargs}` are keyword option-value pairs where the value is a file path the container will need access to. Any keyword that expects a file path should get `--file-` prepended to it. 
    * `{kwargs}` are keyword option-value pairs where the value is anything other than a file path.
    * `--output-keys {keys}` are any output keys as a comma separated list from required the print output statements that contain file paths. `ever_given` automatically handles the file path mounting and unmounting, and will update the file path for us

## Example: Running AutoDock Vina Tutorial Docking Container with `ever_given`
For reference, below is the `--help` info for the Autodock Vina container
```
Options:
  --receptor PATH      path of receptor PDB to dock the ligand into
                       [required]
  --smiles TEXT        SMILES str of ligand to be docked. quote to prevent CLI
                       errors "CCC"  [required]
  --hint PATH          path of hint ligand complex for docking region hint
                       [required]
  --hint_molinfo TEXT  residue name of the ligand in the hint complex
                       [required]
  --hint_radius FLOAT  box size of the box to dock into  [required]
  --output-dir TEXT    Output directory for receptor and docked_ligand files
  --help               Show this message and exit.
```
To run the Autodock Vina container called `adv-tutorial` with `ever_given`, we would change the above options as follows:
   * `{file-kwargs}`: `--hint` and `--receptor` expect file paths. When running with `ever_given`, these would become `--file-hint` and `--file-receptor`.
   * `{kwargs}`: `--hint_radius`, `--hint_molinfo`, and `--smiles` expect values other than file paths. These keyword options are all passed to `ever_given` as is with no changes to the keywords.
   * `--output-keys {keys}`: both Docking output keys `docked_ligand` and `receptor` have file paths as their value. When running with `ever_given` we would add `--output_keys docked_ligand,receptor` to our `python run.py` command.

The final command becomes: `python run.py adv-tutorial --file-hint {pdb_file} --file-receptor {pdb_file} --hint_radius {float} --hint_molinfo {str} --smiles {str} --output-keys docked_ligand,receptor`
