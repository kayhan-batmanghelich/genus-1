import os
import numpy as np
import pandas as import pd
from subprocess import call
import h5py
import argparse
import shutil

def main():

    def snp_overlap(data_dir):
        files_bim = [x for x in os.listdir(data_dir) if \
        '.bim' in x and x[0] != '.']
        chr_list = []
        for i in files_bim:
            try:
                print("Starting {}".format(i))
                chr_list.append(
                    np.genfromtxt(os.path.join(data_dir, i),
                        dtype=str,
                        delimiter='\t',
                        usecols=1))
                print("Finsihed {}".format(i))
            except (IOError, ValueError) as e:
    	        print("{} did not load".format(i))
        intersection = set(chr_list[0]).intersection(*chr_list)
        return list(intersection)

    def reidx(data):
        d = data.copy()
        d.index = [idx for idx in range(data.shape[0])]
        retrn d

    def snp_pgc_match(snp_intersection, pgc_snps, n_snps):
        snp_i = pd.DataFrame(snp_intersection, columns=['snpid'])
        snp_i['snpid_copy'] = snp_i.snpid
        m = pgc_snps.set_index('snpid', 1).loc[snp_i.snpid]
        m = m.dropna(axis=0)
        m['snpid'] = m.index.values
        m = reidx(m)
        matching_psorted = reidx(m.sort_values('p'))
        return matching_psorted.ix[:(n_snps - 1), :]

    def read_bim(data):
        bim = pd.read_csv(data, sep='\t', header=None)
        bim.columns = ['chrom','snpid','gd','pp','a1','a2']
        return bim

    def match(matcher, matchie):
        snpm = matchie.set_index('snpid',1).loc[matcher.snpid]
        snpm['snpid'] = snpm.index.values
        snpm.index= [x for x in range(snpm.shape[0])]
        s = snpm[['chrom', 'snpid', 'gd', 'pp', 'a1', 'a2']]
        return s

    def make_new_bim(snp_reduced, data_dir, out_path):
        files = [x for x in os.listdir(data_dir) \
                if (('bgn' in x) and ('.bim' in x) \
                and (x[0] != '.'))]
        file_store = {fname:None for fname in files}
        for bfile in files:
            try:
                bfile = read_bim(os.path.join(data_dir, bfile))
                match(snp_reduced, bfile).to_csv(
                    '{}'.format(out_path + "/" + bfile),
                    sep='\t', index=None, header=None)
            except IOError:
                pass
        return None

    def move_files(data_dir, new_path):
        files = [x for x in os.listdir(data_dir) if '.bim' not in x]
        for f in files:
            shutil(os.path.join(data_dir, f),
                   os.path.join(new_path, f))
        return None

    def make_hdf(data_dir):
        files = [x[:-4] for x in os.listdir(data_dir)]
        cmd_str = 'python pyPlink.py --bfile {} --out {} --outFormat hdf5'
        for sfile in files:
            cmd = cmd_str.format(
                os.path.join(data_dir, sfile),
                sfile + '_h5')
            call(cmd)
        return None

    def read_fam(data):
        fam = pd.read_csv(data, header=None, sep=' ')
        fam.columns = ['FID', 'IID', 'pid', 'mid', 'sex', 'affected']
        return fam

    def read_hdf(data):
        h5 = h5py.File(data)['genotype']
        h5df = pd.DataFrame(np.transpose(h5))
        return h5df

    def combine_hdf(fam_dir):
        hdf_frames = []
        hdf_dir = os.listdir('.')
        files = [x[:-3] for x in hdf_dir]
        for hdf_file in files:
            hdf = read_hdf(os.path.join(hdf_dr, + hdf_file + '_h5'))
            fam = read_fam(os.path.join(fam_fir, + hdf_file, + '.fam'))
            hdf[id_name] = fam.id_name
            hdf_frames.append(hdf)
        pd.concat(hdf_frames).to_hdf('{}'.format(x[:-3] + '.hdf'))
        return None

    make_new_bim(snp_pgc_match(snp_overlap(snp_path), rsnp, 100000), data_path)
    move_files(data_path, move_path)
    make_hdf(move_path)
    combine_hdf(move_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", help="path to .bim/.bam/.fam files", type=str)
    parser.add_argument("-snp", help="path to snp text file", type=str)
    parser.add_argument("-rsnp", help="path to file to reduce snps", type=str)
    parser.add_argument("-id_name", help="name of ID variable for matching", type=str)
    parser.add_argument("-nbo", help="new bim output path", type=str)
    parser.add_argument("-base_path", help="base path", type=str)
    parser.add_argument("-move_path", help="will create directory to move files to", type=str)
    args = parser.parse_args()
    data_path = args.d
    snp_path = args.snp
    rsnp_path = args.rsnp
    id_name = args.id_name
    nbo = args.nbo
    base_path = args.base_path
    move_path = args.move_path
    main()
