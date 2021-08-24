# Description:
* Docks the ligand represented by the SMILES string into the given receptor using `oechem` and `oedocking`
* Outputs the docked ligand PDB and receptor PDB 

# Setup:
1. `cp oe_license.txt SAMPL-league/attic/oedock-eg`

# Build:
1. `cd SAMPL-league/attic/oedock-base`
2. `docker build -t oedock-base .`
3. `cd SAMPL-league/attic/oedock-eg`
4. `docker build -t oedock-eg .`

# Run:

### Options
```
% docker run -it --rm oedock-eg --help
Usage: oedock [OPTIONS]

Options:
  -s, --smiles TEXT    SMILES string  [required]
  -r, --receptor TEXT  receptor PDB  [required]
  --hint_ligand TEXT   PDB of ligand docked into receptor to improve docking.
                       Must be used without --boxcoords --boxsize and
                       --center.
  --c_x FLOAT          box center x coordinate; must be used with --c_y/z and
                       --sz_x/y/z
  --c_y FLOAT          box center y coordinate; must be used with --c_x/z and
                       --sz_x/y/z
  --c_z FLOAT          box center z coordinate; must be used with --c_x/y and
                       --sz_x/y/z
  --sz_x INTEGER       box size in the x direction; must be used with
                       --c_x/y/z and --sz_y/z
  --sz_y INTEGER       box size in the y direction; must be used with
                       --c_x/y/z and --sz_x/z
  --sz_z INTEGER       box size in the z direction; must be used with
                       --c_x/y/z and --sz_x/y
  --b_xmin FLOAT       box coordinate x min must be used with
                       --b_ymin/zmin/xmax/ymax/zmax
  --b_ymin FLOAT       box coordinate y min must be used with
                       --b_xmin/zmin/xmax/ymax/zmax
  --b_zmin FLOAT       box coordinate z min must be used with
                       --b_xmin/ymin/xmax/ymax/zmax
  --b_xmax FLOAT       box coordinate x max must be used with
                       --b_xmin/ymin/zmin/ymax/zmax
  --b_ymax FLOAT       box coordinate y max must be used with
                       --b_xmin/ymin/zmin/xmax/zmax
  --b_zmax FLOAT       box coordinate z max must be used with
                       --b_xmin/ymin/zmin/xmax/ymax
  --output-dir PATH    Output Directory
  --help               Show this message and exit.
```


### Run commands
`python ../../ever_given/run.py oedock-eg --file-receptor <receptor_file> --smiles <SMILES_str> --sz_x <boxsize_x> --sz_y <boxsize_y> --sz_z <boxsize_z> --c_x <center_x> --c_y <center_y> --c_z <center_z> --output-keys <ligand_key>,<receptor_key>`
* Ex: `python ../../ever_given/run.py oedock-eg --file-receptor ../data/oedock/4w51-cryo.pdb --smiles "CC(C)Cc1ccccc1 " --sz_x 14 --sz_y 14 --sz_z 14 --c_x -32.355 --c_y 7.263 --c_z 2.207 --output-keys docked_ligand,receptor`

`python ../../ever_given/run.py oedock-eg --file-receptor <receptor_file> --smiles <SMILES_str> --b_xmin <boxcoord_xmin> --b_ymin <boxcoord_ymin> --b_zmin <boxcoord_zmin> --b_xmax <boxcoord_xmax> --b_ymax <boxcoord_ymax> --b_zmax <boxcoord_zmax> --output-keys <ligand_key>,<receptor_key>`
* Ex: `python ../../ever_given/run.py oedock-eg --file-receptor ../data/oedock/4w51-cryo.pdb --smiles "CC(C)Cc1ccccc1 "--b_xmin 39.355 --b_ymin 0.2629999999999999  --b_zmin -4.793 --b_xmax -25.354999999999997 --b_ymax 14.263 --b_zmax 9.207 --output-keys docked_ligand,receptor`

`python ../../ever_given/run.py oedock-eg --file-receptor <receptor_file> --smiles <SMILES_str> --file-hint-ligand <hint_ligand_file> --output-keys <ligand_key>,<receptor_key>`
* Ex: `python ../../ever_given/run.py oedock-eg --file-receptor ../data/oedock/4xni_rec.pdb --smiles "CC(C)Cc1ccccc1 " --file-hint-ligand ../data/oedock/4xni_lig.pdb --output-keys docked_ligand,receptor`
