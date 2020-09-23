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
    def __init__(self, target_type, query_type, prody_obj, resnum):
        if not 'COMBS' in os.environ:
            print('Please define COMBS path in os environment.')
            sys.exit()
        date = '20181009'
        self.base = os.environ['COMBS'] 
        self.target_type = target_type
        self.query_type = query_type
        self.resnum = resnum
        self.residue = prody_obj.select('resnum {}'.format(self.resnum))
        self.resname = self.residue.getResnames()[0]
        self.df_path = df_path(self.target_type, self.query_type,
                self.resname, base=self.base)
        self.df = truncate_df(pd.read_pickle(self.df_path))
        score_log_likelihood(self.df)
        self.rotation, self.translation = self.get_transform()
        self.df = apply_transformation(self.rotation, self.translation,
                self.df)

    def get_transform(self):
        ala_path = os.path.join(self.base, 'ideal_ala_' +
                self.query_type + '.pkl')
        ala_df = pd.read_pickle(ala_path)
        align_atoms = rel_coords_dict[self.query_type]

        target_coords = residue_coords_from_prody(self.residue,
                self.resnum, align_atoms)

        ala_xyz = []
        for atom in align_atoms:
            row = ala_df[ala_df['name']==atom]
            assert(row.shape[0]==1)
            xyz = [row['c_x'].iloc[0], row['c_y'].iloc[0], row['c_z'].iloc[0]]
            ala_xyz.append(xyz)
        
        xyz_array = np.array(ala_xyz)
        rotation, translation = get_superimpose_transformation(xyz_array, target_coords)
        return rotation, translation

    # @property
    def unique_vdms(self):
        groups = self.df.groupby(['iFG_count', 'vdM_count',
            'cluster_score'])
        return groups.size().reset_index()[['iFG_count', 'vdM_count',
            'cluster_score']]


def df_path(target_type, query_type, resname,
        base = '/Users/codykrivacic/intelligent_design/combs_database_/'):
    base_extended = os.path.join(base, 'representatives', 'hb_only')
    path = os.path.join(base_extended, target_type)
    date = find_df_dates(base)[path]
    return os.path.join(path, date, query_type,
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


def normalize_dataframe_to_residue(target_coords, label, df,
        base='/Users/codykrivacic/intelligent_design/combs_database_/'):
    rotation, translation = get_transform(target_coords, label, base=base)
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
    target = atoms.select('chain A')
    #target = atoms.select('chain B')

    init('-extrachi_cutoff 0 -ex1 -ex2')
    print('IMPORTING POSE')
    pose = pose_from_file('test_inputs/6dkm.pdb')
    print('SPLITTING POSE')
    pose = pose.split_by_chain(2)
    print('POSE SPLIT')

    #interface_residues = [63, 20, 56, 44, 65]
    interface_residues = [65, 58]
    import time
    for resi in interface_residues:
        dataframes = []
        # df_target = pd.DataFrame()
        print('TRANSFORMING DATABASE FOR RESI {}'.format(resi))
        start_time = time.time()
        combs = Combs('carboxamide', 'SC', target, resi)
        total = time.time() - start_time
        print('Time elapsed: {}'.format(total))

        print('Getting unique iFG/vdM combinations')
        unique = combs.unique_vdms()
        # print(unique)
        # print(unique[(unique['iFG_count']=='1001') &
            # (unique['vdM_count']=='2')])
        print('# unique vdMs: {}'.format(unique.shape[0]))

        print('Iterating through vdMs')
        start_vdm = time.time()
        # df_vdm = pd.DataFrame()
        for index, row in unique.iterrows():
            print('ANALYZING VDM {} OF {}'.format(index, unique.shape[0]))
            iFG = row['iFG_count']
            vdM = row['vdM_count']
            invrot = InvRot(combs, iFG, vdM)
            invrot.create_inverse_rotamers()
            start_resi = time.time()
            for query_residue in range(1, pose.size() + 1):
                df_resi = pd.DataFrame(invrot.get_superimpose_transformations(pose,
                        query_residue))
                #df_vdm = pd.concat([df_vdm, df_resi], ignore_index=True)
                df_resi['iFG_count'] = iFG
                df_resi['vdM_count'] = vdM
                df_resi['cluster_score'] = row['cluster_score']
                df_resi['target_resi'] = resi
                dataframes.append(df_resi)
            end_resi = time.time() - start_resi
            print('Time to iterate all residues: {}'.format(end_resi))
            # df_vdm['iFG_count'] = iFG
            # df_vdm['vdM_count'] = vdM
            # df_vdm['cluster_score'] = row['cluster_score']
            # df_target = pd.concat([df_target, df_vdm], ignore_index=True)
        # df_target['target_resi'] = resi
        df_out = pd.concat(dataframes, ignore_index=True)
        df_out.to_csv(
                os.path.join(
                        'test_outputs', 
                        'df_target_{}.csv'.format(resi)
                        )
                )
        end_vdm = time.time() - start_vdm
        print('Time for vdM iteration: {}'.format(end_vdm))
