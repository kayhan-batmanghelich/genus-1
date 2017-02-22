## scripts:
### match.py
    from match import Match
    # load your data first
    b = pd.read_csv('brain/braindata.csv')
    c = pd.read_csv('cog/cogdatacsv')
    g = pd.read_hdf('gene/genedata.hdf)
    m = Match('IID', [b,c,g])
    brain, cog, gene = m.fit()
    # or
    m = Match('IID', [c,g])
    cog, gene = m.fit()
    
* note that the index will be the ID

### summary.py
    from summary import Summary
    s = Summary(brain, ('group', 'sex'))
    results = s.amount()
    
