import os
import numpy as np
from nipype import Function
from nipype import Node
from nipype import Workflow
from nipype import IdentityInterface

ds="/storage/gablab001/data/genus/fs_cog/pred_diag/data_sets"
data_sets = [os.path.join(ds, x) for x in os.listdir(ds) if ".csv" in x]
response_var = os.path.join(ds, "response.txt")

wf = Workflow(name="classify_disease")
wf.base_dir = "/om/scratch/Sat/ysa"

Iternode = Node(IdentityInterface(fields=['data', 'classifier']), name="Iternode")
Iternode.iterables = [
     ('data', data_sets), 
     ('classifier', ['et', 'lg'])
]

def run(data, classifier, response):
    import numpy as np
    import pandas as pd
    from custom import Mods
    from custom import utils
    
    y = np.genfromtxt(response)
    X = pd.read_csv(data)
    data_mod = data.split('/')[-1].replace('.csv', '')    

    if classifier == 'et':
        od = '/storage/gablab001/data/genus/fs_cog/pred_diag/extra_trees/results/'
        on = classifier + '_{}.pkl'.format(data_mod)
        mod = Mods.FuzzyTrees(X=X, y=y, out_dir=od, out_name=on)
        outpath = mod.run()
    elif classifier == 'lg':
        od = '/storage/gablab001/data/genus/fs_cog/pred_diag/lg/results/'
        on = classifier + '_{}.pkl'.format(data_mod)
        mod = Mods.Logistic(X=X, y=y, out_dir=od, out_name=on)
        outpath = mod.run()

    return outpath

Run = Node(interface=Function(
       input_names = ['data','classifier','response'],
       output_names = ['outpath'],
       function = run), name = 'Run'
)

Run.inputs.response = response_var
wf.connect(Iternode, 'data', Run, 'data')
wf.connect(Iternode, 'classifier', Run, 'classifier')
sbatch_params = '--mem=4G -t 5-23:00:00 --qos=gablab'
wf.run(plugin='SLURM', plugin_args={'sbatch_args': sbatch_params})
 
