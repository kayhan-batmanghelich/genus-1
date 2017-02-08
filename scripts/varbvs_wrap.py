from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu

def run(infile, outfile, colnum):
    import nipype.interfaces.matlab as Matlab
    import os
    def outnames(col, out):
        return out + '{}.mat'.format(col)
    matlab = Matlab.MatlabCommand()
    matlab.inputs.paths = ['varbvs-MATLAB']
    matlab.inputs.script = """
    load('{}', 'I', 'G', 'Z', 'colnames');
    G = single(G);
    yI=I(:,{});
    sigma=(0.2:0.08:1);
    sa=(0.025:0.025:0.4);
    fit=varbvs(G,Z,yI,colnames,[],struct('logodds',-5:0.25:-3));
    w=normalizelogweights(fit.logw);
    PIP=fit.alpha * w(:);
    mu=fit.mu;
    alpha=fit.alpha;
    save('{}','PIP','w','alpha','mu');
    """.format(infile, int(colnum), outnames(colnum, outfile))
    res = matlab.run()

Run = pe.Node(name = 'Run',
	interface = Function(input_names = [
		'infile', 'outfile', 'colnum'],
	output_names = [''],
	function = run))

Iternode = pe.Node(niu.IdentityInterface(fields=['colnum']), name = 'Iternode')

def csv(colnum, outfile):
    import pandas as pd
    import os
    df = pd.DataFrame(columns=['colNum', 'matFn'])
    for i in range(1, colnum):
        df.loc[i] = [i, outfile+'{}.mat'.format(i)]
    df.iloc[:,0] = df.iloc[:,0].astype(int)
    df.to_csv('BFRESULTSFILELIST.csv', index=None)
	
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--infile', type=str, help='input file and path')
parser.add_argument('-o', '--outfile', type=str, help='output path and file')
parser.add_argument('-c', '--cols', type=int, help='number of columns for brain regions')
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
cols = args.cols
csv(cols, outfile)
Run.plugin_args = {'bsub_args': '-q big'}
Run.inputs.outfile = outfile
Run.inputs.infile = infile
Iternode.iterables = [('colnum', [x for x in range(1, cols+1)])]
wf = pe.Workflow(name='updated_bf')
wf.connect(Iternode, 'colnum', Run, 'colnum')
wf.base_dir = '/data/petryshen/yoel/variational_bayes_wrap/'
wf.run(plugin='LSF', plugin_args={'bsub_args': '-q big'})
