import os
import json
import numpy as np
import pandas as pd

def loadjson(x):
    with open(x) as data_file:
        data = json.load(data_file)
    return data

def return_df(sub_id, stat_file, fs_dir, fjson_obj):
    data = np.genfromtxt(
        os.path.join(fs_dir.format(sub_id),
        stat_file), dtype=object)
    col_heads = fjson_obj[stat_file][0].split()
    df = pd.DataFrame(data = data)
    if df.shape[1] == 1:
        df = df.T
    df.columns = col_heads
    return df

def check_data_type(data, column):
    copy = data.copy()
    for idx in range(copy[column].shape[0]):
        try:
            val = int(copy[column][idx])
            copy[column][idx] = val
        except ValueError:
            try:
                val = float(copy[column][idx])
                copy[column][idx] = val
            except ValueError:
                pass
    return copy

def change_df(df):
    for column in df.columns:
        df = check_data_type(df, column)
    return df

headers = loadjson('scripts/statsheaders.json')
os1 = np.genfromtxt('ds000115_R2.0.0/controls.txt', dtype=str)
os2 = np.genfromtxt('ds000115_R2.0.0/schizs.txt', dtype=str)
fs_dir = '/om/user/ysa/open_scz/scz_subjects_dir/{}/stats'

# schiz = 1, control = 0
open_neurosubs = [i for x in (os1, os2) for i in x]
sub_type = lambda x, y: 0 if x in y else 1
opensubs = {i: sub_type(i, os1) for i in open_neurosubs}

huge_data_frame = []
for stat_file in headers.keys():
    print("DOING {}".format(stat_file))
    print("\n")
    single_sub = return_df(opensubs.keys()[0],
    stat_file, fs_dir, headers)
    m = {i:[] for i in single_sub['StructName'].values}

    for sub in opensubs.keys():
        data = change_df(return_df(sub, stat_file, fs_dir, headers))
        data['subid'] = [sub for i in range(data.shape[0])]
        groups = data.groupby('StructName')
        for name, group in groups:
            m[name].append(group)

    for structure in m.keys():
        sdf = pd.concat(m[structure])
        sdf = sdf.reset_index()
        cols = [sdf['StructName'][0] + '_' + i for i in sdf.columns]
        sdf.columns = cols
        sdf = sdf.drop(structure + '_' + 'StructName', 1)
        sdf.to_csv('stats_dataframes/{}.csv'.format(structure), index=None)
        huge_data_frame.append(sdf)

bdf = pd.concat(huge_data_frame, axis=1).to_csv('stats_dataframes/combined_stats.csv', index=None)
