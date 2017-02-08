import pandas as pd
import numpy as np

d_atlas = pd.read_csv('GENUS_FS_ATLAS_D.csv', low_memory=False)
by_study = d_atlas.groupby('STUDY')

def get_group_n(data, group):
    return (data.GROUP == group).sum()

# name, shape, sczcntshape, sczcnt%shape, scz%, cnt%
to_fill = np.zeros((len(by_study), 6), dtype=object)

for idx, (name, group) in enumerate(by_study):
    scz = get_group_n(group, 'Schizophrenia')
    cnt = get_group_n(group, 'Control')
    gN = scz + cnt
    to_fill[idx][0] = name
    to_fill[idx][1] = group.shape[0]
    to_fill[idx][2] = gN
    to_fill[idx][3] = gN / float(group.shape[0])
    to_fill[idx][4] = scz / float(gN)
    to_fill[idx][5] = cnt / float(gN)

pd.DataFrame(
    data=to_fill,
    columns=[
        'Site', 'Total_Samples',
        'SCZ_CNT_Samples', 'Percent_use',
        'SCZ_%', 'CNT_%'
    ]
).to_csv("site_makeup.csv", index=None)

