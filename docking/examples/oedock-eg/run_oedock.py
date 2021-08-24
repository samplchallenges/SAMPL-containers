# imports for python library
import os
import os.path
import shutil
import tempfile

# imports for python packages
from openeye import oechem, oedocking
import openeye.oequacpac as oequacpac
import openeye.oeomega as oeomega
from openeye import oespruce
from openmoltools import openeye
import click

from ligand import Ligand
from receptor import Receptor
from docking import Docking


LIGAND_KEY = "docked_ligand"
RECEPTOR_KEY = "receptor"
SCORE_KEY = "score"



@click.command()
@click.option("-s","--smiles",default="CCC",required=True,help="SMILES string")
@click.option("-r","--receptor",required=True,help="receptor PDB")
@click.option("--hint",help="PDB of ligand docked into receptor to improve docking. Must be used without --boxcoords --boxsize and --center.")
@click.option("--hint_molinfo")
@click.option("--hint_radius")
@click.option("--output-dir")
def oedock(smiles, receptor, hint, hint_molinfo, hint_radius, output_dir) -> None:
	print("beginning", flush=True)
	tempdir = tempfile.mkdtemp() 

	# make a ligand from smiles string and write to oeb
	ligand_oeb_file = os.path.join(tempdir,"ligand_prepped.oeb")
	lig = Ligand(smiles)
	lig.generate_conformers()
	lig.charge()
	lig.write(ligand_oeb_file, oechem.OEFormat_OEB)


	print("making receptor from pdb", flush=True)
	# make receptor from pdb and write to oeb
	boxcoords = None
	receptor_oeb_file = os.path.join(tempdir,"receptor.oeb")
	rec = Receptor(receptor)
	rec.set_make_receptor_param(hint, boxcoords)
	rec.make_receptor()
	rec.write(receptor_oeb_file)



	dock_method = oedocking.OEDockMethod_Chemgauss4
	dock_resolution = oedocking.OESearchResolution_Default

	dock = Docking(dock_method, dock_resolution)
	dock.read_receptor_file(receptor_oeb_file)
	dock.initialize()

	num_poses = 1
	print("docking ligand", flush=True)
	dock.dock(ligand_oeb_file,num_poses)
	print("docked ligand", flush=True)


	docked_lig_pdb_file = f"{output_dir}/docked_ligand.pdb"
	score = dock.write_docked_ligands(docked_lig_pdb_file, oechem.OEFormat_PDB)

	docked_rec_pdb_file = f"{output_dir}/docked_receptor.pdb"
	dock.write_receptor(docked_rec_pdb_file, oechem.OEFormat_PDB)

	print(f"{LIGAND_KEY} {docked_lig_pdb_file}")
	print(f"{RECEPTOR_KEY} {docked_rec_pdb_file}")
	print(f"{SCORE_KEY} {score}")


	shutil.rmtree(tempdir)


if __name__ == "__main__":
	oedock()
