

import argparse
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--chr", "-c", dest="chrom", help="chromosome number", type=str, required=True)
req_grp.add_argument("--effect_i", "-e", dest="effect_i", help="effect iteration", type=str, required=True)

#req_grp.add_argument("--legend", "-l", dest="legend", help="SNP legend file", type=str, required=True)
#req_grp.add_argument("--outpre", "-o", dest="outpre", help="prefix for output file - should be informative about mating", type=str, required=True)
#req_grp.add_argument("--effects", "-e", dest="effects", help="file containing effect sizes", type=str, required=True)
#req_grp.add_argument("--mates", "-m", dest="mates", help="file with mate pair information", type=str, required=True)
#parser.add_argument("--iteration", "-i", dest="iteration", help="iteration", type=str, default="1", nargs="?")
args=parser.parse_args()

import allel
import pandas as pd
import numpy as np
from scipy.stats import binom
import random

random.seed(args.effect_i)

print("reading files")

#load genotypes
ht_sibs = np.load("test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr"
                    +args.chrom
                    +"_e1_sibs_assort.haps.npz")

ht_sibs = allel.HaplotypeArray(ht_sibs['arr_0'])

#load genetic effect sizes simulated in the training dataset
legend = np.loadtxt("test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr"
                      +args.chrom
                      +"_e1_sibs_assort.legend",
                      skiprows=1,dtype="str"
                     )

#load genetic effect sizes simulated in the training dataset
effects = pd.read_csv("train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train."
                        +args.effect_i
                        +".thinned_100kb.effects",
                      header = None,
                      sep = "\t",
                      names = ["ID","A1","effect"]
                     )

#load the .sample file which contains mate pair IDs in the 5th column in the form ("mom,dad")
pairs_np = np.loadtxt("test/sibs/genotypes/genos_complex_assort_1_4.sample",
                        skiprows=1,dtype="str")

var_id = legend[:,0]

print("generating scores")
#to multiply the effects to genotypes, we need both to be in the same order
#generate a list of corresponding row indices for effects data.frame and genotype data
ix_list=[]
for i in range(0,len(effects)):
   if effects.ID.values[i] in var_id:
       ix_list.append((i,np.where(var_id ==effects.ID.values[i])[0][0]))


#subset variants from genotype data and effects file in the correct order
ht_causal = ht_sibs[[item[1] for item in ix_list],:]

#calculate no. of alt alleles for each variant since that's all that's needed for PRS calculation
ht_causal = ht_causal.to_genotypes(ploidy=2)
ht_causal = ht_causal.sum(axis=2)

effects_causal = effects.effect.values[ [item[0] for item in ix_list] ]

#the dot product of the two arrays will yield the PRS for each individual
gsum = effects_causal.dot(ht_causal)

#write prs to file
gsum = np.column_stack( (pairs_np[:,0:2],
                         [args.chrom]*len(gsum),
                         gsum) )

np.savetxt("test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr"
            +args.chrom+"_e"
            + args.effect_i
            +"_sibs_assort.sscore",
           gsum,
           header="FID IID CHROM SCORESUM",
           delimiter=" ",
           fmt="%s",
          comments='')
