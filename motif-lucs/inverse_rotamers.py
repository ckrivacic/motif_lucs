from pyrosetta import *
#import pyrosetta.rosetta
import numpy as np
from utils import *
from numeric import *
from itertools import compress
from database import rel_coords_dict
from contact_labels import *


"""
Generate inverse rotamers based on functional group placement in combs
database, then attach to design chain and get rotation and translation
compared to some 'base' level.
"""


class InvRot(object):

    def __init__(self, combs, ifg, vdm):
        self.combs = combs
        self.ifg = ifg
        self.vdm = vdm
        #self.query = combs.query
        self.alignment_atoms = rel_coords_dict[self.combs.target_type]
        self.vdm_type = 'hb_sc'
        self.get_rows()
        self.vdm_atoms = all_interactions[self.vdm_type][self.vdm_restype]
        self.get_motif_positions()
        #self.pdb_path = 'test_inputs/6dkm.pdb'
        # constraints = 'test_inputs/test_constraints.cst'
        #self.constraints = 'test_inputs/8cho_cst_E.cst'

        #self.pose = rosetta.core.import_pose.get_pdb_and_cleanup(self.pdb_path)

        # Maybe move these outside of class for optimization
        dummy_pose = rosetta.core.pose.Pose()
        self.global_residue_type_set = dummy_pose.residue_type_set_for_pose()
        # self.residue_type = global_residue_type_set.get_representative_type_base_name("ASP")

        #self.restraints = parse_restraints(self.constraints)

    def get_rows(self):
        self.sub_df = self.combs.df[(self.combs.df['iFG_count']==str(self.ifg)) &
                (self.combs.df['vdM_count']==str(self.vdm))] 
        self.vdm_restype = self.sub_df.iloc[0]['resname']

    def get_motif_positions(self):
        sub_df = self.sub_df[self.sub_df['name'].isin(self.vdm_atoms)]
        self.positions = []
        for idx, row in sub_df.iterrows():
            position = {}
            position['atom_name'] = row['name']
            position['coord'] = [row['c_x'], row['c_y'], row['c_z']]
            self.positions.append(position)

    def create_inverse_rotamers(self):
        """
        Creates a set of inverse rotamers that are overlaid onto a set of
        coordinate restraints. 

        In the future, the positioning information coming from the motif
        database will be relative. To address this, I'll need 2
        transformation steps; one that aligns the motif coordinates onto my
        protein, then another to create the inverse rotamers.
        """
        restype = self.global_residue_type_set.get_representative_type_base_name(
                        self.vdm_restype
                )
        print("creating rotamer set")
        rotamer_set = bb_independent_rotamers_extra_chi(restype)

        xyzarray = []
        restrained_atoms = []
        for position in self.positions:
            xyzarray.append(position['coord'])
            restrained_atoms.append(position['atom_name'])

        xyzarray = np.array(xyzarray)

        for residue in rotamer_set:
            residue_restrained_xyz = []
            residue_all_xyz = []
            for atom in restrained_atoms:
                residue_restrained_xyz.append(np.array(residue.xyz(atom)))
            residue_restrained_xyz = np.array(residue_restrained_xyz)
            #for index in range(1, residue.natoms() + 1):
            #    residue_all_xyz.append(np.array(residue.xyz(index)))
            #residue_all_xyz = np.array(residue_all_xyz)

            rotation, translation = \
                    get_superimpose_transformation(residue_restrained_xyz, xyzarray)

            # The following transformation should be saved in rotamer_set,
            # so that's ultimately what we'll return.
            residue.apply_transform_Rx_plus_v(np_array_to_xyzM(rotation),\
                    np_array_to_xyzV(translation))

        self.rotamer_set = rotamer_set

    def get_superimpose_transformations(self, pose, resi):
        '''Get rotation and translation matrices of a pose's residue to
        each inverse rotamer '''
        atoms = ['N', 'C', 'CA']
        rows = []
        # query_resi = 65
        for rot in self.rotamer_set:
            row = {
                    # 'query_resi': query_resi,
                    'query_resi': resi,
                    }
            xyz_rot = []
            xyz_pose = []
            for atom in atoms:
                xyzV = xyzV_to_np_array(rot.xyz(atom))
                xyz_rot.append(xyzV)
                xyzV_pose = xyzV_to_np_array(pose.residue(resi).xyz(atom))
                xyz_pose.append(xyzV_pose)
            xyz_rot = np.array(xyz_rot)
            xyz_pose = np.array(xyz_pose)
            rotation, translation = \
                    get_superimpose_transformation(xyz_pose, xyz_rot)
            row['rot'] = rotation
            row['trans'] = translation
            rows.append(row)

        return rows

#cst_test = ConstrainToInvRot()
#rotamer_set = cst_test.create_inverse_rotamers('GLU')
#cst_test.choose_rotamer()
#cst_test.make_constraints_from_inverse_rotamer()

#designable, repackable = choose_designable_residues(cst_test.pose, [38])
#task_factory = setup_task_factory(cst_test.pose, designable, repackable,
#        motif_dict={38:'E'},layered_design=False, prepare_focus=True)

#bb_movable = [i for i in range(1,cst_test.pose.size() + 1)]
#sc_movable = []
#fast_relax(cst_test.pose, bb_movable, sc_movable, selectors=False)
#print(cst_test.pose.constraint_set())
#print(designable)
#loops = generate_loops_from_res_selector(cst_test.pose, designable, 38)
#fast_design(cst_test.pose, designable, repackable, task_factory=task_factory)


#model_loops(cst_test.pose, designable, repackable, 38,
#        task_factory=task_factory, 
#        fast=True, mover='lhk', resbuffer=4)

#lhk = lhk_xml(task_factory)
#lhk.apply(cst_test.pose)


if __name__=="__main__":
    init("-ignore_unrecognized_res -extrachi_cutoff 0 -ex1 -ex2 -out:overwrite " +\
            "-lh:db_path=/home/ckrivacic/rosetta/database/loophash_db/ " +\
            "-lh:loopsizes 6")
