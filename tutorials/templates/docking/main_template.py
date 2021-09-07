import click
import os.path

@click.command()
@click.option("--receptor", required=True, type=click.Path(exists=True), help="path of receptor PDB to dock the ligand into")
@click.option("--smiles", required=False, type=click.Path(exists=True), help="file with SMILES strings of ligands to be docked")

@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint_molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint_radius",required=True,type=float,help="box size of the box to dock into")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

def docking_main(receptor, smiles,  hint, hint_molinfo, hint_radius, output_dir):
        ''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:        file    receptor PDB path to dock ligand into
                       smiles:          str     SMILES string of ligand to be docked 
                       hint:            file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo:    str     resname of the ligand used in the hint PDB
                       hint_radius:     float   radius around the hint ligand to consider in docking
                       output_dir:      str     output director for receptor and docked_ligand
            OUTPUTS:   prints           docked_ligand {path_to_docked_ligand_file}
                       prints           receptor {path_to_receptor_file}
                       writes file(s)   docked ligand file as a .pdb .mol2 or .sdf
                       writes file(s)   receptor prepped and used by program in docking as .pdb
        '''
        
        # set the output file names / paths
        docked_ligand_file_name = ""
        receptor_file_name = ".pdb"
        path_to_docked_ligand_file = os.path.join(output_dir, docked_ligand_file_name)
        path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
        
	
	
	# YOUR CODE GOES HERE 
        
	
        
        # write out the docked ligand file to path_to_docked_ligand_file
        # write out the receptor file to path_to_receptor_file
       

        # print out the key value pairs 
        #    * where the keys are docked_ligand and receptor
        #    * where the values are the file paths of the docked_ligand and receptor
        
        print(f"docked_ligand {path_to_docked_ligand_file}")
        print(f"receptor {path_to_receptor_file}")


if __name__ == "__main__":
	docking_main()
