from openeye import oechem

def write_OEMol(outstream: str, oformat: "oechem.OEFormat", mol: oechem.OEGraphMol):
	''' opens the oechem.oemolostream() with the path outstream
		
	    Parameters:
	    ----------
	    outstream: str
	               String representing the path of the outstream

	    Returns:
	    --------
	    ofs: oechem.oemolostream 
		 outstream to write the molecules to
	'''
	ofs = oechem.oemolostream()
	ofs.SetFormat(oformat)
	if not ofs.open(outstream):
		oechem.OEThrow.Fatal(f"Unable to open {outstream} for writing")
	oechem.OEWriteMolecule(ofs, mol)


def ReadFromPDB(pdb_file, mol):
	ifs = oechem.oemolistream()
	ifs.SetFlavor(oechem.OEFormat_PDB, oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa

	if not ifs.open(pdb_file):
		oechem.OEThrow.Fatal("Unable to open %s for reading." % pdb_file)

	temp_mol = oechem.OEGraphMol()
	if not oechem.OEReadMolecule(ifs, temp_mol):
		oechem.OEThrow.Fatal("Unable to read molecule from %s." % pdb_file)
	ifs.close()
    
	fact = oechem.OEAltLocationFactory(temp_mol)
	mol.Clear()
	fact.MakePrimaryAltMol(mol)
	return (mol)
