# genus

## For the scripts:
### match.py
    from math import Match
    b = 'brain/GENUS_FS_ATLAS_D.csv'
    g = 'genomic/h5_100k_eur_combined.hdf'
    c = 'cognitive/GENUS_neuropsych_data_Domain_Scores.csv'
    m = Match(b, c, g, 'IID')
    brain, cog, gene = Match.out()
    
