import openeye.oechem as oechem
import openeye.oedocking as oedocking

import oe_io


class Receptor():
	def __init__(self, pdb: str):
		self._bio_system = oechem.OEGraphMol()
		input_success = oe_io.ReadFromPDB(pdb, self._bio_system)
		if not input_success:
			oechem.OEThrow.Fatal("Unable to input protein from PDB file.")
		self._receptor = None
		self._third_make_receptor_param = None


	def make_receptor(self):
		self._receptor = oechem.OEGraphMol()
		oedocking.OEMakeReceptor(self._receptor, self._bio_system, self._third_make_receptor_param)
		


	def set_make_receptor_param(self, hint_ligand, boxcoords):
		make_rec_param = None
		if hint_ligand != None:
			bound_ligand = oechem.OEGraphMol()
			ifs = oechem.oemolistream()
			ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)
			if not ifs.open(hint_ligand):
		                oechem.OEThrow.Fatal("Unable to open %s for reading." % hint_ligand)
			if not oechem.OEReadMolecule(ifs, bound_ligand):
				oechem.OEThrow.Fatal("Unable to open %s for reading." % hint_ligand)
			make_rec_param = bound_ligand

		elif boxcoords:
			xmin,ymin,zmin,xmax,ymax,zmax = boxcoords
			make_rec_param = oedocking.OEBox(xmin, ymin, zmin, xmax, ymax, zmax)

		else:
			oechem.OEThrow.Fatal("Unable to get receptor parameter")

		self._third_make_receptor_param = make_rec_param


	def write(self, receptor_oeb_file: str):
		oedocking.OEWriteReceptorFile(self._receptor, receptor_oeb_file)


	@staticmethod
	def get_boxcoords(boxsize: (int,int,int), center: (float,float,float)):
		center_x, center_y, center_z = center
		size_x, size_y, size_z = boxsize

		xmin = center_x - (0.5 * size_x)
		ymin = center_y - (0.5 * size_y)
		zmin = center_z - (0.5 * size_z)
		
		xmax = center_x + (0.5 * size_x)
		ymax = center_y + (0.5 * size_y)
		zmax = center_z + (0.5 * size_z)

		return (xmin, ymin, zmin, xmax, ymax, zmax)	

