import click
import os.path


# the following are decorators related to the MAIN FUNCTION, these '@click.command()' and '@click.option()' decorators
# MUST remain directly above your main function
@click.command()
@click.option("--receptor", required=True, type=click.Path(exists=True), help="path of receptor PDB to dock the ligand into")
@click.option("--smiles", required=False, type=click.Path(exists=True), help="file with SMILES strings of ligands to be docked")

@click.option("--hint",required=True,type=click.Path(exists=True),help="path of hint ligand complex for docking region hint")
@click.option("--hint_molinfo",required=True,help="residue name of the ligand in the hint complex")
@click.option("--hint_radius",required=True,type=float,help="box size of the box to dock into")

@click.option("--output-dir",help="Output directory for receptor and docked_ligand files")

# @click.option("--your_argument", type=click.Path(exists=True), help="Any special file arguments you")
def docking_main(receptor, smiles,  hint, hint_molinfo, hint_radius, output_dir) -> None:
        ''' docks the given smiles string into the receptor within the area specified by hint and hint-radius
            INPUTS:    receptor:        file    receptor PDB path to dock ligand into
                       smiles:          str     SMILES string of ligand to be docked 
                       hint:            file    hint PDB contains a receptor ligand complex to show binding site region
                       hint_molinfo:    str     resname of the ligand used in the hint PDB
                       hint_radius:     float   radius around the hint ligand to consider in docking
                       output_dir:      str     output director for receptor and docked_ligand\
		       
            OUTPUTS:   prints           docked_ligand {path_to_docked_ligand_file}
	    				where the path_to_docked_ligand_file = output_dir/docked_ligand_file
                       prints           receptor {path_to_receptor_file}
		       			where the path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
                       writes file      docked ligand file as a .pdb .mol2 or .sdf
                       writes file      receptor prepped and used by program in docking as .pdb
		       
	    RETURNS:   None
	    	       Please note that anything your function returns will be ignored by our automated scoring 
		       All outputs (docked_ligand and receptor files) MUST BE saved on disk to the specified
		       output_dir AND the absolute path must be printed out as specified in "OUTPUTS"
        '''
        
        # set the output file names / paths
        docked_ligand_file_name = ""
        receptor_file_name = ".pdb"
        path_to_docked_ligand_file = os.path.join(output_dir, docked_ligand_file_name)
        path_to_receptor_file = os.path.join(output_dir, receptor_file_name)
        
	
	
	# YOUR CODE GOES HERE 
        
	
        
        # write out the docked ligand file to path_to_docked_ligand_file, your docked ligand file must be 
	# saved to the path specified by the output_dir parameter 
	
        # write out the receptor file to path_to_receptor_file, your receptor file must be 
	# saved to the path specified by the output_dir parameter 
       

	# Your final ligand and receptor files should be SAVED to the output_dir (specified as a parameter). 
	# Your main function should also PRINT out the 'key value' pairs 
	#	* key: either the 'docked_ligand' or 'receptor'. 
	# 	* value: the absolute file path to the file on disk. The absolute file paths to your files should be specified
	#	         as 'output_dir/your_file' where output_dir is the output_dir (specified as a parameter). Your final
	# 		 output files MUST BE saved to the output_dir, otherwise our automated scoring will not be able to 
	#		 find them
	# 		* path_to_docked_ligand_file = os.path.join(output_dir,"rec-dock.pdb") 
	# 		* path_to_receptor_file = os.path.join(output_dir,"best_dock.pdb")
        
        print(f"docked_ligand {path_to_docked_ligand_file}")
        print(f"receptor {path_to_receptor_file}")
	
	# A NOTE ABOUT RETURN STATEMENTS
	# Anything your container returns will be ignored. Please make sure that any outputs follow the format
	# in the previous comment


if __name__ == "__main__":
	docking_main()
