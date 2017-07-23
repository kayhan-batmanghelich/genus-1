import os
import numpy as np
from nipype import Function, Node, Workflow, IdentityInterface

wf = Workflow(name='brain_bcv')
wf.base_dir = "/om/user/ysa"

Iternode = Node(IdentityInterface(fields=['col_idx', 'cv_idx']), name = 'Iternode')
Iternode.iterables = [('col_idx', np.arange(170) + 1), ('cv_idx', np.arange(10) + 1)]

def cv_maker(data_path, save_path):
    import scipy.io
    from sklearn.model_selection import StratifiedKFold
    X = scipy.io.loadmat(data_path)['I']
    y = scipy.io.loadmat(data_path)['response'][0]
    cv = StratifiedKFold(n_splits=10, random_state=1)
    train_idx, test_idx = {}, {}
    for idx, (train, test) in enumerate(cv.split(X, y)):
        train_idx['train_{}'.format(idx + 1)] = train + 1
        test_idx['test_{}'.format(idx + 1)] = test + 1
    scipy.io.savemat(save_path, mdict={"train":train_idx, "test":test_idx})
    return save_path
 
CV_maker = Node(interface=Function(
    input_names = ['data_path', 'save_path'],
    output_names = ['save_path'],
    function = cv_maker
), name = 'CV_maker')

CV_maker.inputs.data_path = "/storage/gablab001/data/genus/brain_genomic_bayes/brain_gene.mat"
CV_maker.inputs.save_path = "/storage/gablab001/data/genus/brain_genomic_bayes/cv_idxs.mat"

def run_bayes(in_file, cv_file, cv_idx, col_idx, out_file):
    import cPickle as pickle
    import numpy as np
    import os
    import nipype.interfaces.matlab as Matlab
    def outnames(col, out):
        return os.path.join(out, '{}.mat'.format(col))
    col_names = np.genfromtxt("/storage/gablab001/data/genus/brain_genomic_bayes/170_columns.txt", dtype=str)
    col_save_name = col_names[col_idx - 1] + "_{}_{}_BF".format(cv_idx, col_idx)
    with open("/storage/gablab001/data/genus/brain_genomic_bayes/bayes_reg.m", "r") as src:
        script = src.read().replace("\n", "")
    mat_file = outnames(in_file[:-4] + col_save_name, out_file)
    matlab = Matlab.MatlabCommand()
    matlab.inputs.paths = [
    '/storage/gablab001/data/genus/current/variational_bayes_wrap/varbvs/varbvs-MATLAB',
    '/storage/gablab001/data/genus/current/variational_bayes_wrap/varbvs',
    '/storage/gablab001/data/genus/current/variational_bayes_wrap/varbvs/varbvs-R']
    matlab.inputs.script = script.format(in_file, cv_file, cv_idx, col_idx, mat_file)
    res = matlab.run()
    return mat_file

Run_bayes = Node(interface=Function(
    input_names = ['in_file', 'cv_file', 
                   'cv_idx', 'col_idx',
                   'out_file'],
    output_names = ['mat_file'],
    function = run_bayes
), name='Run_bayes')

Run_bayes.inputs.in_file = "/storage/gablab001/data/genus/brain_genomic_bayes/brain_gene.mat"
Run_bayes.inputs.out_file = "/storage/gablab001/data/genus/brain_genomic_bayes/"

wf.connect(CV_maker, 'save_path', Run_bayes, 'cv_file')
wf.connect(Iternode, 'cv_idx', Run_bayes, 'cv_idx')
wf.connect(Iternode, 'col_idx', Run_bayes, 'col_idx')
wf.run(plugin='SLURM', plugin_args={'sbatch_args':'--mem=12G -t 3-23:00:00', 'max_jobs': 100})
