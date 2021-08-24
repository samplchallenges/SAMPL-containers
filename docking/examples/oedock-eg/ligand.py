import openeye.oechem as oechem
import openeye.oeomega as oeomega
import openeye.oequacpac as oequacpac
import oe_io

class Ligand():
	def __init__(self, smiles: "SMILES str"):
		self._smiles = smiles
		self._oemol = oechem.OEMol()
		oechem.OEParseSmiles(self._oemol, self._smiles)


	def get_smiles(self):
		return self._smiles

	def get_oemol(self):
		return self._oemol

	def write(self, outfile: str, oeformat: "oechem.OEFormat"):
		oe_io.write_OEMol(outfile, oeformat, self._oemol)

	def generate_conformers(self, maxconfs: int=100):
		# Generage conformers
		omega = oeomega.OEOmega()
		omega.SetMaxConfs(100)
		omega.SetStrictStereo(False)
		omega(self._oemol)

	def charge(self, charge_method=oequacpac.OEAM1BCCELF10Charges()):
		oequacpac.OEAssignCharges(self._oemol, oequacpac.OEAM1BCCELF10Charges())

	def reset():
		self._oemol = oechem.OEMol()
		oechem.OEParseSmiles(self._oemol, self._smiles)
