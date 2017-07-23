import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import scipy.io
sns.set_style('whitegrid')

cnames = {
'aliceblue':            '#F0F8FF',
'antiquewhite':         '#FAEBD7',
'aqua':                 '#00FFFF',
'aquamarine':           '#7FFFD4',
'azure':                '#F0FFFF',
'beige':                '#F5F5DC',
'bisque':               '#FFE4C4',
'black':                '#000000',
'blanchedalmond':       '#FFEBCD',
'blue':                 '#0000FF',
'blueviolet':           '#8A2BE2',
}

snpsloc = pd.read_csv('bf_100snp_output/100k_snp_loc_pval.csv')
toplot = scipy.io.loadmat('bf_100snp_output/100snp_170b_32.mat')
snpsloc['PIP'] = toplot['PIP']
snpsloc['-log10PIP'] = np.log10(snpsloc['PIP'])
colors = np.random.choice(cnames.keys(),22)
df = snpsloc.copy()
# 'chromosome','loci','pvalue','snp','alpha'
df.columns = ['chromosome', 'pvalue', 'snpid', 'PIP', 'log10PIP']
df.chromosome = df.chromosome.astype('category')
df = df.sort_values('chromosome')
df['ind'] = range(len(df))
df_grouped = df.groupby(('chromosome'))

def plot_manhattan(data):
    fig = plt.figure(figsize=(20,10))
    fig.suptitle('Bayes Factor Manhattan Plot', fontsize=20)
    ax = fig.add_subplot(111)
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(data):
        group.plot(
            kind = 'scatter', 
            x = 'ind', 
            y = 'log10PIP',
            colors = colors[num % len(colors)],
            ax = ax
        )
        if np.any(group['PIP'] > .5):
            x_labels.append(name)
        else:
            x_labels.append('')
        x_labels_pos.append((
            group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2
        ))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df)])
    ax.set_ylim([-4.4, .5])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('Chromosome', fontsize=20)
    plt.ylabel('log10PIP', fontsize=20)
    plt.savefig('manhattan_plot.svg', dpi=300, bbox_inches='tight')
    plt.close(fig)

plot_manhattan(df_grouped)
