
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--vcf_file","-v",dest="vcf",help="vcf file path",type=str,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file path",type=str,required=True)
args=parser.parse_args()

import numpy as np
import pandas as pd
import allel

vcf = allel.read_vcf(args.vcf)
gt =  allel.GenotypeArray(vcf['calldata/GT'])

pops = [list(range( ((i-1)*250)+1 , ((i-1)*250)+251)) for i in range(1,35)]

a,b,c =allel.weir_cockerham_fst(gt, pops)
a = np.take(a,indices=1,axis=1)
b = np.take(b,indices=1,axis=1)
c = np.take(c,indices=1,axis=1)

denom = a+b+c
ix_0 = np.where(denom!=0)[0]
fst = np.zeros(len(a))
fst[ix_0] =  a[ix_0]/denom[ix_0]

df = pd.DataFrame({"a": a,
                   "b": b,
                   "c": c,
                   "fst": fst})

df.to_csv(args.outpre+".scikit.fst",sep="\t",header=False,index=False)
