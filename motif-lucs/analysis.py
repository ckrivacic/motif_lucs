import pandas as pd
import numpy as np
import sys, os, glob
from ast import literal_eval
# import re

def bin_df(df, degrees=10, angstroms=1):
    '''Bin a dataframe by rotation & translation.'''
    # print(df['rot'])
    # print(type(df['rot'].values[0]))

    nrbins = int(360 // degrees)
    tmax = df['tmax'].max()
    tmin = df['tmin'].min()
    ntbins = int((tmax - tmin) // angstroms)
    rbins = np.linspace(-180, 180, nrbins+1)
    tbins = np.linspace(df['tmin'].min(), df['tmax'].max(), ntbins + 1)
    def bin_rot(array):
        return bin_array(array, rbins)

    def bin_trans(array):
        return bin_array(array, tbins)

    def bin_array(array, bins):
        '''Digitize a numpy array'''
        inds = np.digitize(array, bins)
        binned = tuple([bins[inds[n]-1] for n in range(array.size)])
        return binned

    df['rgroup'] = df['rot'].apply(bin_rot)
    df['tgroup'] = df['trans'].apply(bin_trans)
    
    groups = df.groupby(['rgroup', 'tgroup', 'query_resi', 'target_resi'])

    l = []
    i = 0
    for name, group in groups:
        i +=1
        row = {
                'rot': name[0],
                'trans': name[1],
                'query_resi': name[2],
                'target_resi': name[3],
                'size': len(group.index)
                }
        l.append(row)
        #print(name)
        #print(group)
    l = pd.DataFrame(l).sort_values('size', ascending=False,
            ignore_index=True)
    print(l.iloc[0:50])
    print(l[l['rot']==(0,0,0)])
    print(i)
    l.to_pickle('binned.pkl')


def prep_df(path):
    # df = pd.read_csv(sys.argv[1], converters={'rot': literal_eval}, nrows=10000)
    from scipy.spatial.transform import Rotation as R
    def make_rot_array(text):
        text = text.replace('[', '')
        text = text.replace(']', '')
        text = text.replace('\n', ' ')
        array = np.fromstring(text, sep=' ')
        array = array.reshape(3, 3)
        array = R.from_matrix(array)
        array = array.as_euler('xyz', degrees=True)
        return array

    def make_trans_array(text):
        text = text.replace('[', '')
        text = text.replace(']', '')
        text = text.replace('\n', ' ')
        array = np.fromstring(text, sep=' ')
        return array

    if path.endswith('csv'):
        df = pd.read_csv(path, nrows=10000)
        df['rot'] = df['rot'].apply(make_rot_array)
        df['trans'] = df['trans'].apply(make_trans_array)
    elif path.endswith('pkl'):
        df = pd.read_pickle(path)
    df['tmin'] = df['trans'].apply(lambda x: min(x))
    df['tmax'] = df['trans'].apply(lambda x: max(x))

    return df

def prep_folder(path):
    dfs = []
    for filename in glob.glob('{}/*.pkl'.format(path)):
        dfs.append(prep_df(filename))
    df = pd.concat(dfs, ignore_index=True)

    return df

if __name__=='__main__':
    path = sys.argv[1]
    print('Loading dataframe(s)')
    if os.path.isfile(path):
        df = prep_df(sys.argv[1])
    else:
        df = prep_folder(path)
    print('Dataframes loaded; binning')
    bin_df(df, degrees=30, angstroms=10)
