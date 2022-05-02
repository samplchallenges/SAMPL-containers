from rdkit import Chem
from rdkit.Chem import AllChem
import os
import os.path
import mdtraj as md
import tempfile
import shutil
from collections import namedtuple



DockBox = namedtuple("DockBox", ["c_x", "c_y", "c_z", "sz_x", "sz_y", "sz_z"])

class Autodock():
    def __init__(self, python_path: str, utility_path: str, vina_path: str, receptor: str, hint:str=None, molinfo:str=None, radius:int=None):
        self.python_path = python_path         
        self.utility_path = utility_path
        self.vina_path = vina_path

        if hint != None and os.path.exists(hint) and radius != None and radius >= 0 and molinfo != None:
            box = self._get_hint_box(hint, molinfo, radius)

        self.receptor_file_path = receptor
        assert(os.path.exists(self.receptor_file_path))
        self.dockbox = box


    def __str__(self):
        return f"Autodock: Receptor=\"{self.receptor_file_path}\", Box=({self.sz_x},{self.sz_y},{self.sz_z},{self.c_x},{self.c_y},{self.c_z}), Python=\"{self.python_path}\", Utilities=\"{self.utility_path}\", Vina=\"{self.vina_path}\""


    def _get_hint_box(self, hint_pdb: str, molinfo: str, radius: int):
        ''' get the docking box so Autodock knows where to dock the ligand
            we calculate the dockbox by finding the center of mass from the 
            hint ligand and building a box in 1 radius above, below, and side
            side to side from the the center of mass
        '''
        temp_dir = tempfile.mkdtemp()
        hint_lig_pdb = os.path.join(temp_dir,"hint_lig.pdb")
        Autodock.extract_ligand(hint_pdb, molinfo, hint_lig_pdb)
        hint_lig = md.load(hint_lig_pdb)
        com = md.compute_center_of_mass(hint_lig)
        c_x = com[0][0]*10
        c_y = com[0][1]*10
        c_z = com[0][2]*10
        sz_x = 2*radius
        sz_y = 2*radius
        sz_z = 2*radius

        shutil.rmtree(temp_dir)
        return DockBox(c_x, c_y, c_z, sz_x, sz_y, sz_z)
        
    def _get_utility_cmd(self, pyfile: str):
        ''' returns the python run command for the given pyfile
        '''
        pyfile_path = os.path.join(self.utility_path, pyfile)
        return f"{self.python_path} {pyfile_path}"

    def _make_config_file(self, exhaustiveness, num_modes, flex, config_path):
        ''' makes the configuration file which is an input to the vina 
            command
        '''
        configfile = open(config_path, "w")
        configfile.write(f"receptor = {self.recprep_pdbqt}\n")

        configfile.write(f"center_x = {self.dockbox.c_x}\n")
        configfile.write(f"center_y = {self.dockbox.c_y}\n")
        configfile.write(f"center_z = {self.dockbox.c_z}\n")

        configfile.write(f"size_x = {self.dockbox.sz_x}\n")
        configfile.write(f"size_y = {self.dockbox.sz_y}\n")
        configfile.write(f"size_z = {self.dockbox.sz_z}\n")

        if exhaustiveness != None:
            configfile.write(f"exhaustiveness = {exhaustiveness}\n")
        if num_modes != None:
            configfile.write(f"num_modes = {num_modes}")
        
        configfile.close()


    def _save_highest_score(self, ligdock_path: str, lighighscore_path: str):
        ''' writes the highest scoring docked pose into its own file
        '''

        infile = open(ligdock_path, 'r')
        outfile = open(lighighscore_path, 'w')
        
        for line in infile.readlines():
            if "MODEL" in line:
                continue
            if "ENDMDL" in line:
                break
            outfile.write(line)
                
        outfile.close()
        infile.close()

    def _get_score(self, score_path):
        ''' returns the score output by autodock
        '''
        with open(score_path) as scoref:
            for line in scoref:
                if "Aff" in line:
                    return float(line.split()[1])

    @staticmethod
    def charge_ligand(smiles, outfile: str):
        ''' creates a molecule from the smiles string, charges the ligand
            and adds 3D coordinates
        '''
        mol = Chem.MolFromSmiles(smiles)
        mol2 = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol2)

        # save the charged molecule as an sdf
        ostream = Chem.SDWriter(outfile)
        ostream.write(mol2)

    @staticmethod
    def extract_ligand(complex_pdb, lig_name, out_pdb):
        ''' extract a ligand from a PDB and produce a new PDB with just
            the ligand
        ''' 
        if not os.path.exists(complex_pdb):
            return False
        outfile = open(out_pdb, "w")
        with open(complex_pdb, "r") as complex_file:
            for line in complex_file:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:21].strip() == lig_name:
                    outfile.write(line)
            outfile.write("END")
        outfile.close()
        return True


    @staticmethod
    def sdf_to_pdbqt(sdf_file: str, pdbqt_file: str):
        ''' static method to convert sdf file to pdbqt file using
            openbabel
        '''
        assert(os.path.exists(sdf_file))
        os.system(f"obabel {sdf_file} -O {pdbqt_file}")

    @staticmethod
    def pdbqt_to_pdb_static(pdbqt_path, pdb_path):
        ''' static method to convert pdbqt to pdb file using cut
        '''
        os.system(f"cut -c-66 {pdbqt_path} > {pdb_path}")

    def pdbqt_to_pdb(self, pdbqt_path, pdb_path):
        ''' convert pdbqt file to pdb file using MGLTools
        '''
        cmd = self._get_utility_cmd("pdbqt_to_pdb.py")
        os.system(f"{cmd} -f {pdbqt_path} -o {pdb_path}")

    def prep_ligand(self, ligchg_pdbqt: str, ligprep_pdbqt: str):
        ''' prepare ligand PDBQT for autodock
        '''
        cmd = self._get_utility_cmd("prepare_ligand4.py")
        os.system(f"{cmd} -l {ligchg_pdbqt} -o {ligprep_pdbqt}")

    def prep_receptor(self, recprep_pdbqt: str):
        ''' prepare receptor PDBQT for autodock
        '''
        cmd = self._get_utility_cmd("prepare_receptor4.py")
        print(cmd)
        self.recprep_pdbqt = recprep_pdbqt
        os.system(f"{cmd} -r {self.receptor_file_path} -o {recprep_pdbqt}")


    def dock(self, ligprep_pdbqt: str, ligdock_pdbqt: str, lighighscore_pdbqt: str, lighighscore_pdb: str, exhaustiveness: int, num_modes: int, flex):
        ''' dock the ligand in the receptor using autodock, return the affinity score 
        '''
        temp_dir = tempfile.mkdtemp()
        config_path = os.path.join(temp_dir,"config.txt")
        score_path = os.path.join(temp_dir,"score.txt")

        self._make_config_file(exhaustiveness, num_modes, flex, config_path)
        os.system(f"{self.vina_path} --cpu 1   --config {config_path} --ligand {ligprep_pdbqt} --out {ligdock_pdbqt}")

        self._save_highest_score(ligdock_pdbqt, lighighscore_pdbqt)
        self.pdbqt_to_pdb(lighighscore_pdbqt, lighighscore_pdb)

        os.system(f"{self.vina_path} --cpu 1 --config {config_path} --ligand {lighighscore_pdbqt} --out /tmp/rescore.pdbqt --score_only | grep 'Affinity:' > {score_path}")
        score = self._get_score(score_path)

        shutil.rmtree(temp_dir)
        return score

    def save_receptor_pdb(recprep_pdb: str):
        ''' save the prepped PDB used for ligand docking, should be in the same frame of 
            any docked ligand PDBs output
        '''
        Autodock.pdbqt_to_pdb(recprep_pdbqt, recprep_pdb)
