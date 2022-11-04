#python built in
import os
import os.path
import random
import sys
import tempfile
import shutil
import csv
import time

# python packages
import click
from rdkit import Chem
from rdkit.Chem import AllChem
from autodock import Autodock,DockBox

SCORE_KEY = "docking_score"
BINDS_KEY = "compound_binds"

PYTHON_PATH = "/opt/app/dependencies/mgl/bin/python"
UTILITIES_PATH = "/opt/app/dependencies/mgl/MGLToolsPckgs/AutoDockTools/Utilities24"
VINA_PATH = "/opt/app/dependencies/adv/bin/vina"

SCORE_CUTOFF = -8.8

def print_debug(debug: bool, msg:str):
	print(f"{msg}\n" if debug else "", end="")
	sys.stdout.flush()

# set command line options
@click.command()
@click.option("--receptor",required=True,type=click.Path(exists=True),help="path of receptor PDB to dock the ligand into")
@click.option("--smiles",required=True,help="SMILES str of ligand to be docked. quote to prevent CLI errors \"CCC\"")

@click.option("--c_x")
@click.option("--c_y")
@click.option("--c_z")
@click.option("--sz_x")
@click.option("--sz_y")
@click.option("--sz_z")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

@click.option('--batch', is_flag=True, help="Input and output will be batch format")

@click.option('--debug', is_flag=True,help="prints debug print statements when --debug flag is used")
def main_function(receptor, smiles, c_x, c_y, c_z, sz_x, sz_y, sz_z, output_dir, batch, debug) -> None:
	''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
			INPUTS:	receptor:	 file	receptor PDB path to dock ligand into
					   smiles:	   str	 SMILES string of ligand to be docked, use quotes
					   hint:		 file	hint PDB contains a receptor ligand complex to show binding site region
					   hint_molinfo: str	 resname of the ligand used in the hint PDB
					   hint_radius:  float   radius around the hint ligand to consider in docking
					   output_dir:   str	 output director for receptor and docked_ligand
					   batch		 bool	bool whether to use regular or batched docking
					   debug:		bool	bool used for degbug print statemetns
	'''
	if not os.path.exists(receptor):
		raise FileNotFoundError(f"File \"{receptor}\" does not exist.")

	dock_box = DockBox(c_x, c_y, c_z, sz_x, sz_y, sz_z)

	if batch:
		docking_batch(receptor, smiles, dock_box, output_dir, debug)

	else:
		docking_no_batch(receptor, smiles, dock_box, output_dir, debug)




def docking_no_batch(receptor, smiles, dock_box, output_dir, debug) -> None:

	# create a temporary directory to store intermediate files
	temp_dir = tempfile.mkdtemp()
	# Set file names for intermediate and output files
	ligchg_sdf_path = os.path.join(temp_dir, "lig-chg.sdf")
	ligchg_pdbqt_path = os.path.join(temp_dir,"lig-chg.pdbqt")
	ligprep_pdbqt_path = os.path.join(temp_dir, "lig-prep.pdbqt")

	receptor_path = receptor
	receptorprep_pdbqt_path = os.path.join(temp_dir, "rec-prep.pdbqt")
	receptor_pdb_path = os.path.join(output_dir,"dock_rec.pdb")

	ligdock_pdbqt_path = os.path.join(temp_dir, "lig-dock.pdbqt")
	highscore_pdbqt_path = os.path.join(temp_dir, "best-dock.pdbqt")
	highscore_pdb_path = os.path.join(output_dir,"dock_best_pose.pdb")

	score = 0
	binds = -1

	try:

		# create Autodock calculator object
		adv = Autodock(PYTHON_PATH, UTILITIES_PATH, VINA_PATH, receptor_path, dock_box=dock_box)

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

		binds = 1 if score < -5 else 0
	except Exception as e:
		print(f"Error occured during docking: {e}")
	finally: 
		print(f"{SCORE_KEY} {score}")
		print(f"{BINDS_KEY} {binds}")

	# Anything your container returns will be ignored. Please make sure that any outputs follow the format
	# in the previous comment


def docking_batch(receptor, smiles, dock_box, output_dir, debug):

	temp_dir = tempfile.mkdtemp()

	ligchg_sdf_path = os.path.join(temp_dir, "lig-chg.sdf")
	ligchg_pdbqt_path = os.path.join(temp_dir,"lig-chg.pdbqt")
	ligprep_pdbqt_path = os.path.join(temp_dir, "lig-prep.pdbqt")

	receptor_path = receptor
	receptorprep_pdbqt_path = os.path.join(temp_dir, "rec-prep.pdbqt")
	receptor_pdb_path = os.path.join(output_dir,"dock_rec.pdb")

	ligdock_pdbqt_path = os.path.join(temp_dir, "lig-dock.pdbqt")
	highscore_pdbqt_path = os.path.join(temp_dir, "best-dock.pdbqt")
	highscore_pdb_path = os.path.join(temp_dir,"dock_best_pose.pdb")

	score_path = os.path.join(temp_dir, "score.txt")

	final_molfile_path = os.path.join(output_dir, "docked_ligands.mol")

	adv = Autodock(PYTHON_PATH, UTILITIES_PATH, VINA_PATH, receptor_path, dock_box=dock_box)
	adv.prep_receptor(receptorprep_pdbqt_path)
	assert(os.path.exists(receptorprep_pdbqt_path))

	Autodock.pdbqt_to_pdb_static(receptorprep_pdbqt_path, receptor_pdb_path)
	fields = ['id', 'name', 'value']

	scores_file = os.path.join(output_dir, "scores.csv")
	scores_f = open(scores_file, "w+")
	scores_writer = csv.writer(scores_f)
	scores_writer.writerow(fields)


	binds_file = os.path.join(output_dir, "binds.csv")
	binds_f = open(binds_file, "w+")
	binds_writer = csv.writer(binds_f)
	binds_writer.writerow(fields)

	with open(smiles, "r") as smiles_fp:
		smiles_reader = csv.DictReader(smiles_fp)
		for row in smiles_reader:
			score = 0
			binds = -1
			try:
				print(f"dock start: {time.time()}")
				molid = int(row['id'])
				name = row['name']
				smiles = row['value']

				print(f"smiles: {smiles}")

				Autodock.charge_ligand(smiles, ligchg_sdf_path)
				assert(os.path.exists(ligchg_sdf_path))
				
				Autodock.sdf_to_pdbqt(ligchg_sdf_path, ligchg_pdbqt_path)
				assert(os.path.exists(ligchg_pdbqt_path))
				adv.prep_ligand(ligchg_pdbqt_path, ligprep_pdbqt_path)

				num_modes = 1
				exhaustiveness = None
				flex = None
				score = adv.dock(ligprep_pdbqt_path, ligdock_pdbqt_path, highscore_pdbqt_path, highscore_pdb_path, exhaustiveness, num_modes, flex)
				print(f"dock end: {time.time()}")

				print(f"score: {score}")
				binds = 1 if score < SCORE_CUTOFF else 0
			except Exception as e:
				print(f"Error occurred during docking: {e}")
			finally:
				scores_writer.writerow([molid, name, score])
				binds_writer.writerow([molid, name, binds])


	shutil.rmtree(temp_dir)

	binds_f.close()
	scores_f.close()

	os.system(f'cat {binds_file}')
	os.system(f'cat {scores_file}')

	print(f"{SCORE_KEY} {scores_file}")
	print(f"{BINDS_KEY} {binds_file}")

if __name__ == "__main__":
	main_function() 
