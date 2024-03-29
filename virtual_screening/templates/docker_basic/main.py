#python built in
import os
import tempfile
import shutil

#pip installable
import click

# Global variables specific to the virtual screening challenge
# DO NOT ALTER
SCORE_KEY = "docking_score"
BINDS_KEY = "compound_binds"
BINDER = 1
NON_BINDER = 0

# set command line options
@click.command()
@click.option("--receptor",required=True,type=click.Path(exists=True),help="path of receptor PDB to dock the ligand into")
@click.option("--smiles",required=True,help="SMILES str of ligand to be docked. quote to prevent CLI errors \"CCC\"")

@click.option("--c_x", help="center of docking box x coordinate")
@click.option("--c_y", help="center of docking box y coordinate")
@click.option("--c_z", help="center of docking box z coordinate")
@click.option("--sz_x", help="size of docking box in x direction")
@click.option("--sz_y", help="size of docking box in y direction")
@click.option("--sz_z", help="size of docking box in z direction")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

@click.option('--batch', is_flag=True, help="Input and output will be batch format")
def main_function(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z, output_dir, batch) -> None:
    ''' scores the given smiles string docked into the receptor within the area specified 
        by c_x/y/z, sz_x/y/z
          INPUTS:    
            receptor:     file      receptor PDB path to dock ligand into
            smiles:       str/file  SMILES string of ligand(s) to be docked. str when 1 mol, csv file when many mols
            c_x:          float     center of docking box x coordinate
            c_y:          float     center of docking box y coordinate
            c_z:          float     center of docking box z coordinate
            sz_x:         float     size of docking box in x direction
            sz_y:         float     size of docking box in y direction
            sz_s:         float     size of docking box in z direction
            output_dir:   str       directory to store output files
            batch:        bool      bool whether to use regular or batched docking
    '''
    if not os.path.exists(receptor):
        raise FileNotFoundError(f"File \"{receptor}\" does not exist.")

    if batch:
        run_screening_batch(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z)

    else:
        run_screening_no_batch(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z, output_dir)


def run_screening_no_batch(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z) -> None:
    ''' Some challenges may be run by passing your container a batch of inputs. 
        In this case, we will be evaluating whether many molecules passed as a 
        csv file of smiles strings are binders to the given receptor. 
          INPUTS: 
            receptor, c_x, c_y, c_z, sz_x, sz_y, sz_z (See main_function)
            smiles:  file  csv file in the format "id,name,value" where id=molid, 
                           name=molname, and value=smiles string
          OUTPUTS:
              prints       docking_score {numerical docking score}
              prints       compound_binds {1 (BINDER) if binds; else 0 (NON_BINDER)}

          RETURNS:
            None     you are welcome to have your function return something, but
                     it will not be used for scoring or record keeping purposes

    '''

    # Your docking/scoring functionality goes here

    print(f"{SCORE_KEY} {score}")
    print(f"{BINDS_KEY} {binds}")

    # Anything your container returns will be ignored. Please make sure that any outputs follow the format
    # in the previous comment


def run_screening_batch(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z, output_dir):
    ''' Some challenges may be run by passing your container a batch of inputs. 
        In this case, we will be evaluating whether many molecules passed as a 
        csv file of smiles strings are binders to the given receptor. 
          INPUTS: 
            receptor, c_x, c_y, c_z, sz_x, sz_y, sz_z, output_dir (See main_function)
            smiles:  file  csv file in the format "id,name,value" where id=molid, 
                           name=molname, and value=smiles string
          OUTPUTS:
              prints       docking_score {path_to_docking_scores_csv}
              prints       compound_binds {path_to_binders_csv}
              writes file  scores csv file with columns id, name, value where 
                           id is the molid, name is the molname, and value is 
                           the numerical docking score
              writes file  binders csv file with columns id, name, value where 
                           id is the molid, name is the molname, and value is 
                           the either 1 (BINDER) if the molecule is a binder or
                           0 (NON_BINDER) if the molecule is a non-binder
          RETURNS:
              None     you are welcome to have your function return something, but
                       it will not be used for scoring or record keeping purposes

    '''
    temp_dir = tempfile.mkdtemp()

    # open and prepare the results csv files to write score and binder/non-binder
    # results into
    fields = ['id','name', 'value']
    scores_file = os.path.join(output_dir, "scores.csv")
    scores_f = open(scores_file, "w+")
    scores_writer = csv.writer(scores_f)
    scores_writer.writerow(fields)

    binds_file = os.path.join(output_dir, "binds.csv")
    binds_f = open(binds_file, "w+")
    binds_writer = csv.writer(binds_f)
    binds_writer.writerow(fields)

    # Loop through each molecule smiles strings in the smiles csv file
    #     * score and evaluate each molecule as a BINDER=1 or NON_BINDER=0
    with open(smiles, "r") as smiles_fp:
        smiles_reader = csv.DictReader(smiles_fp)
        for row in smiles_reader:
            # unpack each row in the csv
            molid = int(row['id'])
            name = row['name']
            smiles = row['value']

            bind_condition = ""
            binds = BINDER if bind_condition else NON_BINDER


            # Your docking/scoring functionality goes here


            scores_writer.writerow([molid, name, score])
            binds_writer.writerow([molid, name, binds])

    shutil.rmtree(temp_dir)

    binds_f.close()
    scores_f.close()

    # Print out your results
    # Results should be the results key such as the value SCORE_KEY, docking_score
    # followed by a space and the path to your results file
    print(f"{SCORE_KEY} {scores_file}")
    print(f"{BINDS_KEY} {binds_file}")

