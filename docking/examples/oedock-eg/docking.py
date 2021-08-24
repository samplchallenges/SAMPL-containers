import openeye.oechem as oechem
import openeye.oedocking as oedocking

class Docking():
	def __init__(self, dock_method: "oedocking.OEDockMethod", dock_resolution: "oedocking.OESearchResolution"):
		self._dock_method = dock_method
		self._sdtag = oedocking.OEDockMethodGetName( self._dock_method )
		self._dock_resolution = dock_resolution
		self._oedock = oedocking.OEDock( self._dock_method, self._dock_resolution )

		self._receptor_file = None
		self._receptor = oechem.OEGraphMol()

		self._docked_ligands = list()


	def read_receptor_file(self, receptor_file: str):
		if not oedocking.OEReadReceptorFile( self._receptor, receptor_file ):
			# raise an exception if the receptor file cannot be read
			raise Exception("Unable to read receptor from {0}".format( receptor_file ))


	def initialize(self):
		if not self._oedock.Initialize(self._receptor):
			# raise an exception if the receptor cannot be initialized
			raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))


	def dock(self, ligand_oeb_input: str, num_poses:int=2):
		ifs = oechem.oemolistream( ligand_oeb_input )

		for lig in ifs.GetOEMols():
			dockedLig, score = self._dock_molecule( num_poses, lig )

			self._docked_ligands.append( (dockedLig, score) )


	def write_docked_ligands(self, ligand_outfile: str, oeformat: "oechem.OEFormat"):
		ofs_lig_pdb = oechem.oemolostream( ligand_outfile )
		ofs_lig_pdb.SetFormat(oeformat)

		dockedLig, score = self._docked_ligands[0]
		oechem.OEWriteMolecule(ofs_lig_pdb, dockedLig)

		return score


	def write_receptor(self, receptor_outfile: str, oeformat: "oechem.OEFormat"):
		ofs_rec_pdb = oechem.oemolostream( receptor_outfile )
		ofs_rec_pdb.SetFormat(oeformat)
		oechem.OEWriteMolecule(ofs_rec_pdb, self._receptor)

	def _dock_molecule( self, num_poses: int, mcmol ) -> tuple:
		''' Docks the multiconfomer molecule, with the given number of poses
	            Returns a tuple of the docked molecule (dockedMol) and its score
	            i.e. ( dockedMol, score )
		'''
		dockedMol = oechem.OEMol()

		#Dock the molecule into a given number of poses
		res = self._oedock.DockMultiConformerMolecule(dockedMol, mcmol, num_poses)

		if res == oedocking.OEDockingReturnCode_Success:
			#Annotate the molecule with the score and SDTag that contains the docking method
			oedocking.OESetSDScore(dockedMol, self._oedock, self._sdtag)
			self._oedock.AnnotatePose(dockedMol)
			score = self._oedock.ScoreLigand(dockedMol)
			oechem.OESetSDData(dockedMol, self._sdtag, "{}".format(score))
			return dockedMol, score

		else:
			# raise an exception if the docking is not successful
			raise Exception("Unable to dock ligand {0} to receptor".format( dockedMol ))











