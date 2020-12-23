
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--burden","-b",dest="burden",help="burden array - npz",type=str,required=True)
req_grp.add_argument("--pop","-p",dest="popfile",help="pop file",type=str,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output prefix",type=str,required=True)

args=parser.parse_args()


import numpy as np
import pandas as pd

print("reading files")
#load genotypes
burden = np.load(args.burden)
burden = burden['arr_0']
#load longitude,latitude information
pop = pd.read_csv(args.popfile,sep="\t")

seed=burden[0,0]

print("calculating gini curves")
for i in range(0,100):
    key1='burden_'+str(i)
    pop[key1]=burden[i,2:]

##Calculate no. of demes the variant appears in
gini_b=np.empty((100,36))
bsum=[]
for k in range(0,100):
    key1='burden_'+str(k)

    x=pop.groupby(['deme','longitude','latitude'],as_index=False)[key1].sum()
    x=x[key1]
    bsum.append(np.sum(x))

    x2=np.cumsum(np.sort(x))/np.sum(x)
    gini_b[k,:]=x2

final_output=np.column_stack((np.repeat(seed,100),
                              np.arange(100),
                              bsum,
                              gini_b))

print("writing results")
#output dispersion
np.savetxt(args.outpre,
           final_output,
           delimiter=",",
           fmt='%3.3g')
