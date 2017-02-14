
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 21:48:50 2014

@author: kayhan
"""

import os
import sys
import h5py


import argparse
from argparse import RawTextHelpFormatter

from plinkio import plinkfile

# importing the packages
import  matplotlib.pylab as plt
import  numpy as np

# pandas to open CVS file
import pandas


class genotypeData:
    def __init__(self,bedFn=None, verbose=0):
        if not(bedFn==None):
            self.readPLINKfile(bedFn)
            self.verbose = verbose
            
    def readPLINKfile(self,bedFn):  
        plink_file = plinkfile.open( bedFn )
        if not plink_file.one_locus_per_row( ):
            print( "This script requires that snps are rows and samples columns." )
            exit( 1 )

        sample_list = plink_file.get_samples( )
        locus_list = plink_file.get_loci( )          
        
        # make sure sample and locus
        self.numSubjects = len(sample_list)
        self.numLoci = len(locus_list)

        self.Data = np.zeros( ( len(locus_list) , len(sample_list)  ) , dtype=np.int8)
        cnt = 0
        for locus, row in zip( locus_list, plink_file ):
            self.Data[cnt,:] = np.array(row)
            cnt = cnt + 1
 
        self.snpInfo = pandas.read_csv(bedFn+'.bim',names=['chrNum','snpID','genDist','physLoc','A1','A2'],sep='\t')
        self.subjInfo = pandas.read_csv(bedFn+'.fam',names=['familyID','individualID','paternalID','maternalID','sex','phenotype'],sep=' ')         
         
    def saveGenotype(self, outFile, outFormat):
        if (outFormat=='hdf5'):
            f = h5py.File(outFile, "w")
            dset = f.create_dataset("genotype", (self.numLoci, self.numSubjects), dtype=np.int8)
            dset[...] = self.Data
            f.close()
        else:
            print('this format is not implemented yet !!!')

     
    def summaryInfo(self): 
        # some sanity check
        AACount = sum(self.Data==0)
        AaCount = sum(self.Data==1)
        aaCount = sum(self.Data==2)
        naCount = sum(self.Data==3)
        print("Number of subjects in the dataset is : ", self.numSubjects )
        print("Number of loci in the dataset  : ", self.numLoci)
        print("Number of AA in the datset : ", AACount)
        print("Number of Aa in the datset : ", AaCount)
        print("Number of aa in the datset : ", aaCount)
        print("Number of NA in the datset : ", naCount)
        print 'Here is the ratios   :  ', '#(Aa) / #(AA) :', (0.0 + AaCount)/AACount, '  #(aa) / #(AA) :', (0.0 + aaCount)/AACount, '#NA /#(AA):', (0.0 + naCount)/AACount
        print("Here is the head of snp info : ")
        print(self.snpInfo.head(5))
        print("Here is the head of subject info : ")
        print(self.subjInfo.head(5))

        # visualize part of matrix
        if (self.verbose>0):
            plt.imshow(self.Data[0:400,0:400])
            plt.colorbar()
            plt.draw()
            plt.show()








if __name__ == '__main__':
    #parse input arguments 
    parser = argparse.ArgumentParser(description="""This program reads the plink format ands write all or part of it to a different format. 
                                                    Here is a list of available options : 
                                                       * HDF5 : 
                                                        
                                                    This is an example of how to use it:
                                                    [Name of Method]   --bfile  ADNI_lastQCb37_Chr19   --out  ADNI_lastQCb37_Chr19.h5   --outFormat hdf5  
                                                    
                                                    ** NOTE ** When you read data from hdf format, remeber that the order in PYTHON and MATLAB is different.""", formatter_class=RawTextHelpFormatter )
                                                       
    parser.add_argument('--bfile', help='Specify .bed, .fam and  ', required=True)
    parser.add_argument('--out', help='Specify output filename', required=True)
    parser.add_argument('--outFormat', help='Specify outputformat', required=True,  type=str, choices=['hdf5', 'mat'])
    #parser.add_argument('-s','--subjListFn',help='Text file containing a list of subjects', required=True)
    #parser.add_argument('-i','--inputCSVFile', help='Input CSV file containing the histogram data', required=True)
    #parser.add_argument('-o','--outputCSVFile', help='Output CSV file containing the new features', required=True)
    #parser.add_argument('-O','--outputModelFile', help='Output model file containing the topic model', required=True)
    #parser.add_argument('-v','--verboseLevel', help='Level of verbose: 0: none (default), 1: text, 2: text+graph', required=False, default=0 , type=np.int16 )
    #parser.add_argument('--preProcessSteps', metavar='preProcessSteps', type=str, nargs='+', help='Steps of pre-processing to apply ', choices=['cleanup', 'l2_normalize'])
    args = vars(parser.parse_args())
    print args 
    bedFile = args['bfile']
    outFile = args['out']
    outFormat = args['outFormat']
    
    g = genotypeData(bedFile)
    
    if (outFormat in ['hdf5']):
        g.saveGenotype(outFile,outFormat)     
    else:
        assert 0, "This format is not supported : " + outFormat  
            
    #preProcessSteps = args['preProcessSteps']    
    #inputCSVFile = args['inputCSVFile']
    #outputCSVFile = args['outputCSVFile']
    #verboseLevel = args['verboseLevel']
    #preProcessSteps = args['preProcessSteps']
    #maxIteration = args['maxIteration']
    #outputModelFile = args['outputModelFile']
