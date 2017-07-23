import os
import pandas as pd

sites = [
'scz_camh_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_cati_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_cida_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_cogs_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_deco_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_gapl_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_imh4_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_imhs_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_lanr_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_mcic_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_mfs1_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_mts1_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_nefs_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_paco_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_page_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_phrs_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_tcdn_all_gb-qc.hg19.ch.fl.bgn.reid',
'scz_umcu_all_gb-qc.hg19.ch.fl.bgn.reid',
]

atlas = pd.read_csv('/data/petryshen/yoel/atlas/GENUS_FS_ATLAS_D.csv', low_memory=False)
datadir = '/data/petryshen/yoel/20160823/all/'

def getFam(directory, site):
    f = pd.read_csv(os.path.join(directory, site+'.fam'), header=None, sep=' ')
    f.columns = ['FID', 'IID', 'PaternalID', 'MaternalID', 'sex', 'affected']
    print("{}".format(site))
    return f

def removallNaN(fam, phen, id1, id2):
    try:
        matched = phen.set_index(id1, drop=False).loc[fam[id2]]
        matched.drop_duplicates(inplace=True)
        matched.drop('IID', 1, inplace=True)
        droplist = []
        for idx in matched.index:
            droplist.append(all(pd.isnull(matched.loc[idx])))
        matched_nan_dropped = matched.drop(matched.index[droplist])
        return matched_nan_dropped.shape[0]
    except KeyError:
        print("Error on {} and {}".format(id1, id2))
        return float('NaN')

combinations = [('IID', 'IID'), ('IID', 'FID'), ('FID', 'FID'), ('FID', 'IID')]
combinations_df = pd.DataFrame(index=sites, columns=combinations)

for combo in combinations:
    list_for_df = []
    for i in sites:
        list_for_df.append(removallNaN(getFam(datadir, i), atlas, combo[0], combo[1]))
    combinations_df[combo] = list_for_df

combinations_df.to_csv('fam_with_phen_match.csv')
