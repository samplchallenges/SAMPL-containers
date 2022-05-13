#python built in
import os
import os.path
import random
import sys
import tempfile
import shutil

# python packages
import click
from rdkit import Chem
from rdkit.Chem import AllChem
from autodock import Autodock

LIGAND_KEY = "docked_ligand"
RECEPTOR_KEY = "receptor"

PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"
UTILITIES_PATH = "/opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24"
VINA_PATH = "/opt/app/dependencies/adv/bin/vina"

def print_debug(debug: bool, msg:str):
	print(f"{msg}\n" if debug else "", end="")
	sys.stdout.flush()

# set command line options
@click.command()
@click.option("--receptor",required=True,type=click.Path(exists=True),help="path of receptor PDB to dock the ligand into")
@click.option("--smiles",required=True,help="SMILES str of ligand to be docked. quote to prevent CLI errors \"CCC\"")

@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint-molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint-radius",required=True,type=float,help="box size of the box to dock into")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

@click.option('--debug', is_flag=True,help="prints debug print statements when --debug flag is used")
def main_function(receptor, smiles, hint, hint_molinfo, hint_radius, output_dir, debug) -> None:
	''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:     file    receptor PDB path to dock ligand into
                       smiles:       str     SMILES string of ligand to be docked, use quotes 
                       hint:         file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo: str     resname of the ligand used in the hint PDB
                       hint_radius:  float   radius around the hint ligand to consider in docking
                       output_dir:   str     output director for receptor and docked_ligand
                       debug:        bool    bool used for degbug print statemetns
	'''

	# check to ensure the inputs are valid
	if hint_radius < 0:
		raise ValueError("Hint radius must be a positive number")
	if not os.path.exists(receptor) or not os.path.exists(hint):
		problem_file = receptor if not os.path.exists(receptor) else hint
		raise FileNotFoundError(f"File \"{problem_file}\" does not exist.")

	# create a temporary directory to store intermediate files
	temp_dir = tempfile.mkdtemp()
	# Set file names for intermediate and output files
	ligchg_sdf_path = os.path.join(temp_dir, "lig-chg.sdf")
	ligchg_pdbqt_path = os.path.join(temp_dir,"lig-chg.pdbqt")
	ligprep_pdbqt_path = os.path.join(temp_dir, "lig-prep.pdbqt")

	receptor_path = receptor
	receptorprep_pdbqt_path = os.path.join(temp_dir, "rec-prep.pdbqt")
	receptor_pdb_path = os.path.join(output_dir,"rec-dock.pdb") 

	ligdock_pdbqt_path = os.path.join(temp_dir, "lig_dock.pdbqt")
	highscore_pdbqt_path = os.path.join(temp_dir, "best_dock.pdbqt")
	highscore_pdb_path = os.path.join(output_dir,"best_dock.pdb")

	score_path = os.path.join(temp_dir, "score.txt")

	# create Autodock calculator object
	adv = Autodock(PYTHON_PATH, UTILITIES_PATH, VINA_PATH, receptor_path, hint, hint_molinfo, hint_radius)

	# prepare receptor for docking
	adv.prep_receptor(receptorprep_pdbqt_path)
	assert(os.path.exists(receptorprep_pdbqt_path))
	#assert(os.path.exists(receptorprep_pdbqt_path))


	# convert SMILES str to 3D molecule and prep for docking
	Autodock.charge_ligand(smiles, ligchg_sdf_path)
	assert(os.path.exists(ligchg_sdf_path))

	Autodock.sdf_to_pdbqt(ligchg_sdf_path, ligchg_pdbqt_path)
	assert(os.path.exists(ligchg_pdbqt_path))
	adv.prep_ligand(ligchg_pdbqt_path, ligprep_pdbqt_path)

	# dock the ligand
	num_modes = 1
	exhaustiveness = None
	flex = None
	score = adv.dock(ligprep_pdbqt_path, ligdock_pdbqt_path, highscore_pdbqt_path, highscore_pdb_path, exhaustiveness, num_modes, flex)

	Autodock.pdbqt_to_pdb_static(receptorprep_pdbqt_path, receptor_pdb_path)

	# clean up the temporary directory and intermediate files created 
	shutil.rmtree(temp_dir)

	# Your final ligand and receptor files should be SAVED to the output_dir (specified as a parameter). 
	# Your main function should also PRINT out the 'key value' pairs 
	#	* key: either the LIGAND_KEY 'docked_ligand' or the RECEPTOR_KEY 'receptor'. 
	# 	* value: the absolute file path to the file on disk. The absolute file paths to your files should be specified
	#	         as 'output_dir/your_file' where output_dir is the output_dir (specified as a parameter). Your final
	# 		 output files MUST BE saved to the output_dir
	# 		* receptor_pdb_path = os.path.join(output_dir,"rec-dock.pdb") 
	# 		* highscore_pdb_path = os.path.join(output_dir,"best_dock.pdb")
	print(f"{LIGAND_KEY} {highscore_pdb_path}")
	print(f"{RECEPTOR_KEY} {receptor_pdb_path}")
	
	# Anything your container returns will be ignored. Please make sure that any outputs follow the format
	# in the previous comment

if __name__ == "__main__":
	main_function()
