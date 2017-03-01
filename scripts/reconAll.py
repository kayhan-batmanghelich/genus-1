from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as niu

def run_recon(sub_id):
    import os
    from subprocess import call

    cobre_dir = '/om/user/ysa/open_scz/COBRE'
    open_dir = '/om/user/ysa/open_scz/ds000115_R2.0.0'
    open_path = 'anat/{}_T1w.nii.gz'
    cobre_path = 'session_1/anat_1/mprage.nii.gz'
    sd = '/om/user/ysa/open_scz/scz_subjects_dir'

    if 'sub-' in sub_id:
        t1 = os.path.join(open_dir, sub_id, open_path.format(sub_id))
    else:
        t1 = os.path.join(cobre_dir, sub_id, cobre_path)

    cmd = 'recon-all -i {} -s {} -sd {} -all'.format(t1, sub_id, sd)
    call(cmd.split())
    return None

RunRecon = pe.Node(name = 'RunRecon',
    interface = Function(input_names = ['sub_id'],
                         output_names = [''],
                         function = run_recon))

cobre_dir = '/om/user/ysa/open_scz/COBRE'
open_dir = '/om/user/ysa/open_scz/ds000115_R2.0.0'

cobre_subs = os.listdir(cobre_dir)
open_subs = [x for x in os.listdir(open_dir) if 'sub' in x]

Iternode = pe.Node(niu.IdentityInterface(fields=['sub']), name='Iternode')
Iternode.iterables = [('sub', [i for x in (cobre_subs, open_subs) for i in x])]

wf = pe.Workflow(name='scz_recons')
wf.connect(Iternode, 'sub', RunRecon, 'sub_id')
wf.base_dir = '/om/user/ysa/open_scz'
wf.run(plugin='SLURM', plugin_args={'sbatch_args': '-c2 --mem=12G'})
