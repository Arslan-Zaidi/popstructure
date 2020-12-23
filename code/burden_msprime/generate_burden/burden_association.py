
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--burden","-b",dest="burden",help="burden array - npz",type=str,required=True)
req_grp.add_argument("--pheno","-p",dest="pheno",help="phenotype file",type=str,required=True)
req_grp.add_argument("--cm_pca","-c",dest="cm",help="common PCA file",type=str,required=True)
req_grp.add_argument("--re_pca","-r",dest="re",help="rare PCA file",type=str,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output_prefix",type=str,required=True)

args=parser.parse_args()


import numpy as np
import statistics
import pandas as pd
import statsmodels.api as sm
from scipy import (stats,ndimage)

print("loading files")
#load burden
burden = np.load(args.burden)
burden = burden['arr_0']

#load phenotype file
pheno=pd.read_csv(args.pheno,sep="\t")

#load genetic pcs
gpc_com=np.loadtxt(args.cm,
                       usecols=range(2,102),skiprows=1)

gpc_rare=np.loadtxt(args.re,
                        usecols=range(2,102),skiprows=1)

def assoc(bmat_i,response,gpc_mat,npcs=0):
    #bmat is the burden for the ith gene
    #response is the phenotype (smooth or sharp)
    #npcs is number of PCs to correct for (default=0)
    if(np.isnan(bmat_i).any()):
        return [np.nan,np.nan]
    else:
        if npcs==0:
            predictors=sm.add_constant(bmat_i)
            model=sm.OLS(response,predictors).fit()
            beta=model.params[1]
            pval=model.pvalues[1]

        if npcs>0:
            predictors=np.column_stack((bmat_i,gpc_mat[:,0:(npcs-1)]))
            predictors=sm.add_constant(predictors)
            model=sm.OLS(response,predictors).fit()
            beta=model.params[1]
            pval=model.pvalues[1]

        return [beta,pval]

print("running association between burden and phenotype")
#carry out association test between burden and sharp effect
ngenes=100
np_p1=np.empty((ngenes,8))
seed=burden[0,0]
for i in range(0,ngenes):
        #without correction
        stat_result1=assoc(burden[i,2:],pheno['smooth'],gpc_com,npcs=0)
        stat_result2=assoc(burden[i,2:],pheno['sharp'],gpc_com,npcs=0)

        #correction using 100 common pcs
        stat_result3=assoc(burden[i,2:],pheno['smooth'],gpc_com,npcs=100)
        stat_result4=assoc(burden[i,2:],pheno['sharp'],gpc_com,npcs=100)

        #correction using 100 rare pcs
        stat_result5=assoc(burden[i,2:],pheno['smooth'],gpc_rare,npcs=100)
        stat_result6=assoc(burden[i,2:],pheno['sharp'],gpc_rare,npcs=100)

        np_p1[i,0]=seed
        np_p1[i,1]=i
        np_p1[i,2]=stat_result1[1]
        np_p1[i,3]=stat_result2[1]
        np_p1[i,4]=stat_result3[1]
        np_p1[i,5]=stat_result4[1]
        np_p1[i,6]=stat_result5[1]
        np_p1[i,7]=stat_result6[1]

print("final output, hang tight")

#output pvalue for association test
np.savetxt(args.outpre,np_p1,delimiter=",",fmt='%6.6g')
