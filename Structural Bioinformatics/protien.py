from Bio.PDB import *
import nglview as nv
import ipywidgets
import urllib.request
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils.ProtParam import ProtParamData
from statistics import mean
import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import PPBuilder
import pandas as pd
from treelib import Node, Tree
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 200)


class Quality:
    '''Check Quality of Structure:
    - Resolution of the X-ray data
    - Refinement (the R-factor & Rfree)
    - The B-factors (temperature factors)
    - Ocupency of Atom
    - Model geometry (Ramachandran plot)
    '''
    all = []

    def __init__(self, structure, path):  # innitial value and make it as constructer
        self.structure = structure
        self.path = path
        # assert will be a good way to check what happing matching your expectation

        Quality.all.append(self)

    def poly_detail(self):
        ''' Use this function to build and show polypeptide chain'''
        polypeptide_builder = CaPPBuilder()
        num_of_sequence = 1
        seq = []
        for polypeptide in polypeptide_builder.build_peptides(structure, aa_only=False):
            sequence = polypeptide.get_sequence()
            print(f"Sequence: {num_of_sequence}, Length: {len(sequence)}")
            print(sequence)
            seq.append(sequence)
            num_of_sequence += 1

    def get_resolution(self):
        ''' Get resolution'''
        return structure.header["resolution"]

    def get_r_factor(self):
        ''' Get r factor'''
        try:
            build_dict = MMCIF2Dict(path)  # make a whole mmcif in to dictionary and can use keys to get data
        except ValueError as e:
            print("No R-free value report")
        return build_dict['_refine.ls_R_factor_R_free']

    def get_data(self):
        ''' Get Overall Data'''
        polypeptide_builder = CaPPBuilder()
        poly_list = []
        for index, polypeptide in enumerate(polypeptide_builder.build_peptides(structure, aa_only=False)):
            poly_list.append((index, polypeptide))
        badr = []
        for i in range(0, len(poly_list)):
            for res, (phi, psi) in zip(poly_list[i][1], poly_list[i][1].get_phi_psi_list()):
                badr.append((i + 1, res.id, res.get_resname(), [atom for atom in res],
                             mean([atom.get_bfactor() for atom in res]), mean([atom.get_occupancy() for atom in res]),
                             phi, psi, str(res.get_parent())))

        self.df = pd.DataFrame(badr,
                               columns=['Polychain_number', 'res_id', 'Residue', 'Atoms', 'B-factor', 'Ocupency', 'phi',
                                        'psi', 'Chain_ID'])
        self.b_fac = self.df['B-factor'].mean()
        return self.df

    def cal_b_factor(self):
        ''' Get B-factor by average a whole protien'''
        return self.b_fac

    def cal_bad_b_factor_region(self):
        ''' Show region that high B-factor'''
        mark = self.df['B-factor'] > 80
        x = self.df[mark]
        x['Remark'] = "Bad B-factor region"
        return x[['Polychain_number', 'res_id', 'Residue', 'B-factor', 'Chain_ID','Remark']]

    def ocupency(self):
        ''' Get get ocupency if any atom shows lower ocupency'''
        try:
            mark1 = self.df['Ocupency'] < 1
            y = self.df[mark1]
            y['Remark'] = "low Ocupency"
        except ValueError as e:
            print("No bad occupancy")
        return y[['Polychain_number', 'res_id', 'Residue', 'Ocupency', 'Chain_ID','Remark']]

    def rachamandran_outliers(self):
        ''' Calculate % of Rachamandran disallow region'''
        mask2 = self.df['phi'] > 0
        mask3 = self.df['phi'] < 180
        mask4 = self.df['psi'] < 0
        r = self.df[mask2 & mask3 & mask4]
        return len(r) / len(self.df) * 100

    def rachamandran_outliers_region(self):
        ''' Get data of disallow region'''
        mask2 = self.df['phi'] > 0
        mask3 = self.df['phi'] < 180
        mask4 = self.df['psi'] < 0
        o = self.df[mask2 & mask3 & mask4]
        o['Remark'] = "Disallow Rachmandran"
        return o[['Polychain_number', 'res_id', 'Residue', 'phi', 'psi', 'Chain_ID','Remark']]
    def low_quality_region(self):
        ''' Show low quality region combine from bad b-factor, lower ocupency and
        disallow rachamandran'''
        all_stars = pd.concat([self.cal_bad_b_factor_region(), self.ocupency(),self.rachamandran_outliers_region()])
        return all_stars[['Chain_ID','Polychain_number', 'res_id', 'Residue','Remark']]
####################################################
#  Quality of Score High = 5, Medium = 3, low = 2
####################################################
    def give_score(self):
        ''' Give Quality Score'''
        score = []
        # score from resolution
        if self.get_resolution() < 2.2:
            score.append(5)
        elif 2.2 < self.get_resolution() <= 2.8:
            score.append(3)
        else:
            score.append(1)

        # score from Refinement (the R-factor & Rfree)
        if float(self.get_r_factor()[0]) < 0.22:
            score.append(5)
        elif 0.22 < float(self.get_r_factor()[0]) < 0.4:
            score.append(3)
        else:
            score.append(1)

        # The B-factors (temperature factors)
        if len(self.cal_bad_b_factor_region()) / len(self.get_data()) == 0:
            score.append(5)
        elif len(self.cal_bad_b_factor_region()) / len(self.get_data()) < 20:
            score.append(3)
        else:
            score.append(1)
        # Ocupency of Atom
        if len(self.ocupency()) == 0:
            score.append(5)
        else:
            score.append(3)
        # Model geometry (Ramachandran %of residue in disallow region)
        if self.rachamandran_outliers() < 10:
            score.append(5)
        elif 10 < self.rachamandran_outliers() < 15:
            score.append(3)
        else:
            score.append(1)
        sum(score)
        judge = ""
        if sum(score) > 21:
            judge = "High"
        elif sum(score) > 17 and sum(score) < 21:
            judge = "Medium"
        elif sum(score) < 17:
            judge = "Low"
        return (sum(score),judge)



if __name__ == "__main__":

    # Automatically download structure
    # Download protien file from PBD
    pdbl = PDBList()
    pdbl.retrieve_pdb_file('2WXW')
    parser = MMCIFParser()
    path = 'wx/2WXW.cif'
    structure = parser.get_structure('PHA-L', path)  # Pars structure
    Q = Quality(structure, path)

    print("Poly-peptide detials\n")
    Q.poly_detail()
    print("\n")
    print("Protien structure resolution: ", Q.get_resolution(), "\n")
    print("R-free: ", Q.get_r_factor(), "\n")
    Q.get_data()
    print("Overall B-Factor:", Q.cal_b_factor(), "\n")

    print(f"with bad B-factor residue:")
    print(Q.cal_bad_b_factor_region(), "\n")

    print(f"Ocupency data:")
    print(Q.ocupency(), "\n")

    print("Rachamandran outlier:")
    print(Q.rachamandran_outliers(), "\n")

    print("Ressidues in disallow region:")
    print(Q.rachamandran_outliers_region())


    print("\n")
    print("Total Quality Score: ",Q.give_score())
    print("\n")
    print("Polypeptide regions that are potentially low quality: ")
    print(Q.low_quality_region())

    print("Compare Quality : High to Low ")
    score_lise = [("1FAT", 19), ("3LZT", 23), ("2WXW", 15)]
    score_lise.sort(key=lambda a: a[1], reverse=True)
    tree = Tree()
    tree.create_node(score_lise[0][0], score_lise[0][0])  # root
    for i in range(1, len(score_lise)):
        if i == len(score_lise) - 1:
            tree.create_node(score_lise[i][0], score_lise[i][0], parent=score_lise[i - 1][0])
        elif score_lise[i][1] > 20:
            tree.create_node(score_lise[i][0], score_lise[i][0], parent=score_lise[0][0])
        elif 20 >= score_lise[i][1] > score_lise[len(score_lise) - 1][1]:
            tree.create_node(score_lise[i][0], score_lise[i][0], parent=score_lise[0][0])
    tree.show()