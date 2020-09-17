'''H-bonded sidechain interactions'''
hb_sc = {
        'ARG': ['NH1', 'NH2', 'NE'],
        'GLN': ['NE2', 'CD', 'OE1'],
        'ASN': ['ND2', 'CG', 'OD1'],
        'LYS': ['NZ', 'CE', 'CD'],
        'SER': ['OG', 'CB', 'CA'],
        'TYR': ['OH', 'CZ', 'CE2'],
        'TRP': ['CD1', 'NE1', 'CE2'],
        'HIS': ['NE2', 'CE1', 'ND1'],
        'ASP': ['OD1', 'CG', 'OD2'],
        'GLU': ['OE1', 'CD', 'OE2'],
        }


all_interactions = {
        'hb_sc': hb_sc
        }
