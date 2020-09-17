'''
Functions for interacting with the combs database. 
Here are some things we'll need:
    - Take a functional group and find a match (possibly based on
      Boltzmann-weighted probability distribution)
    - Align database to particular functional group or residue
    - Also, pre-process dataframe. For each residue on target interface,
      save a dataframe with translated/rotated coordinates. Dataframe
      should only cover that type of interaction.
'''

import pandas as pd
import numpy as np
from inverse_rotamers import *
from align import get_superimpose_transformation
import os, prody


rel_coords_dict = dict(SC=['CA', 'N', 'C'], HNCA=['N', 'H', 'CA'],
                       CO=['C', 'O', 'CA'], PHI_PSI=['CA', 'N', 'C'])


class Combs(object):
    def __init__(self, query_type, target_type, prody_obj, resnum):
        date = '20181009'
        self.base = '/Users/codykrivacic/intelligent_design/combs_database_/' 
        self.query_type = query_type
        self.target_type = target_type
        self.resnum = resnum
        self.residue = prody_obj.select('resnum {}'.format(self.resnum))
        self.resname = self.residue.getResnames()[0]
        self.df_path = df_path(self.query_type, self.target_type,
                self.resname)
        self.df = truncate_df(pd.read_pickle(self.df_path))
        score_log_likelihood(self.df)

    def get_transform(self):
        ala_path = os.path.join(self.base, 'ideal_ala_' +
                self.query_type + '.pkl')
        ala_df = pd.read_pickle(ala_path)
        align_atoms = rel_coords_dict[self.query_type]

        query_coords = residue_coords_from_prody(self.residue,
                self.resnum, align_atoms)

        ala_xyz = []
        for atom in align_atoms:
            row = ala_df[ala_df['name']==atom]
            assert(row.shape[0]==1)
            xyz = [row['c_x'].iloc[0], row['c_y'].iloc[0], row['c_z'].iloc[0]]
            ala_xyz.append(xyz)
        
        xyz_array = np.array(xyz_array)
        rotation, translation = get_superimpose_transformation(xyz_array, query_coords)
        return rotation, translation


def df_path(query_type, target_type, resname,
        base = '/Users/codykrivacic/intelligent_design/combs_database_/'):
    base_extended = os.path.join(base, 'representatives', 'hb_only')
    path = os.path.join(base_extended, query_type)
    date = find_df_dates(base)[path]
    return os.path.join(path, date, target_type,
                resname.upper() + '.pkl')


def find_df_dates(base='/Users/codykrivacic/intelligent_design/combs_database_'):
    import glob
    path = os.path.join(base, 'representatives', 'hb_only')
    folders = glob.glob(path + '/*')
    date_dict = {}
    for folder in folders:
        date = os.listdir(folder)[0]
        date_dict[folder] = date
    return date_dict


def truncate_df(df):
    df = df[['iFG_count', 'vdM_count', 'resname', 'resname_vdm', 'name',
        'c_x', 'c_y', 'c_z', 'label', 'cluster_number', 'cluster_size',
        'cluster_type']]
    return df

def score_log_likelihood(df):
    groups = df.groupby(['iFG_count', 'vdM_count'])
    sizes = []
    for name, group in groups:
        sizes.append(group['cluster_size'].iloc[0])
    mean = sum(sizes) / len(sizes)

    log_likelihood = -np.log(df['cluster_size'] / mean)
    df['cluster_score'] = log_likelihood


def select_plane(vdm_num, ifg_num, label, df):
    sub_df = df[(df['iFG_count']==ifg_num) & (df['vDM_count']==vdm_num)]


def normalize_dataframe_to_residue(query_coords, label, df,
        base='/Users/codykrivacic/intelligent_design/combs_database_/'):
    rotation, translation = get_transform(query_coords, label, base=base)
    df = apply_transformation(rotation, translation, df)




def apply_transformation(rotation, translation, df):
    '''Maybe don't do this for the dataframe, and instead transform pose
    onto the ideal alanine? Except then it would be hard to bin. Let's
    just see how long this will take.'''
    def transform_row(row):
        xyz = np.array([row['c_x'], row['c_y'], row['c_z']])
        new_xyz = np.dot(rotation, xyz) + translation
        row['c_x'] = new_xyz[0]
        row['c_y'] = new_xyz[1]
        row['c_z'] = new_xyz[2]
        return row

    df = df.apply(transform_row, axis=1)
    return df


def residue_coords_from_prody(atoms_obj, residue, atoms, chain=None):
    '''Gets coordinates from a prody object for a given residue. Can
    specify chain but best practice is probably to do that once and pass
    only the chain you want to this function since this function will
    probably be called many times.'''
    if chain:
        atoms_obj = atoms_obj.select('chain {}'.format(chain))
    resi = atoms_obj.select('resnum {}'.format(residue))

    xyz = []
    for atom in atoms:
        xyz.extend(resi.select('name {}'.format(atom)).getCoords())

    # Should I make xyz a numpy array before returning?
    return xyz


if __name__=='__main__':
    #combs = Combs('carboxamide', 'SC', 'asn')
    print('Imported dataframe')

    ex_target_path = 'test_inputs/6dkm.pdb'
    atoms = prody.parsePDB(ex_target_path)
    query = atoms.select('chain A')
    #target = atoms.select('chain B')

    interface_residues = [63, 20, 56, 44, 65]
    for resi in interface_residues:
        #xyz = residue_coords_from_prody(query, resi, )
        combs = Combs('carboxamide', 'SC', query, resi)
    xyz = []
    for atom in rel_coords_dict['SC']:
        xyz.extend(tar_sel.select('name {}'.format(atom)).getCoords())
    print('Coordinates for transformation:')
    coords = xyz
    print(coords)

    rot, trans = get_transform(coords, 'SC')
    import time
    start_time = time.time()
    print(combs.df.iloc[1])
    combs.df = apply_transformation(rot, trans, combs.df)
    print(combs.df.iloc[1])
    total = time.time() - start_time
    print('Time elapsed: {}'.format(total))

    init('-extrachi_cutoff 0 -ex1 -ex2')
    invrot = InvRot(combs, 3702, 1)
    invrot.create_inverse_rotamers()
    print('IMPORTING POSE')
    pose = pose_from_file('test_inputs/6dkm.pdb')
    print('SPLITTING POSE')
    pose = pose.split_by_chain(1)
    print('POSE SPLIT')
    start_time = time.time()
    rows = invrot.get_superimpose_transformations(pose, 5)
    print(rows)
    total = time.time() - start_time
    print('Time elapsed: {}'.format(total))
