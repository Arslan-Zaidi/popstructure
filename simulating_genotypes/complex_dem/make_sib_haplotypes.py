


import argparse
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--vcf", "-v", dest="vcf", help="vcf file. needed only for its header", type=str, required=True)
req_grp.add_argument("--mates", "-m", dest="mates", help="file with mate pair information", type=str, required=True)
req_grp.add_argument("--outpre", "-o", dest="outpre", help="prefix for output file - should be informative about mating", type=str, required=True)
parser.add_argument("--iteration", "-i", dest="iteration", help="iteration", type=str, default="1", nargs="?")
args=parser.parse_args()

import allel
import pandas as pd
import numpy as np
from scipy.stats import binom
import random


#load vcf file
vcf = allel.read_vcf(args.vcf)

#get genotypes from vcf
gt = allel.GenotypeArray(vcf['calldata/GT'])

#convert genotypes to haplotypes
#this is necess# Make mate pairs
ht=gt.to_haplotypes()

#get sample ID from vcf
samples = vcf["samples"]

nsnps = len(vcf['variants/POS']) # no. of variants
nsamples = len(samples) # no. of samples


# Now, loop over the pairs and construct haplotypes for siblings by sampling one haplotype from mom and one from dad at random
# instead of doing this with actual haplotypes, which takes up time and memory, I'll do this using just the haplotype indices
# the haplotype indices will then be used to slice columns from ht in that order

#load the .sample file which contains mate pair IDs in the 3rd column in the form ("mom,dad")
pairs = np.loadtxt(args.mates,skiprows=1,dtype="str")

#keep the first half - the pairs are the same for the 2nd sibling
npairs = len(samples)//2
pairs = pairs[0:npairs,]

#separate mom and dad's ID
pairs = np.array([item.split(",") for item in pairs[:,2]])

#get the index of each parent in the samples array
pairs_ix = []
for i in range(0,npairs):
    mom_ix = np.where(samples == pairs[i][0])[0][0]
    dad_ix = np.where(samples == pairs[i][1])[0][0]
    pairs_ix.append((mom_ix,dad_ix))

#loop over each pair and generate two sibling haplotypes
sibs_ix = np.empty((npairs,4),dtype="i4")   #empty array to store sibs indices in

for i in range(0, npairs):
    parent_pair = pairs_ix[i]

    # mom & dad's position in the genotype array (range: 0-8999)
    mom_a = parent_pair[0]
    dad_a = parent_pair[1]

    # mom & dad's haplotype index (0 or 1)
    mom_b = np.random.choice([0,1], size = 2, replace=True)
    dad_b = np.random.choice([0,1], size = 2, replace=True)

    # from the above, generate mom & dad's position in the haplotype array (range: 0-17999)
    mom_ht = 2*mom_a + mom_b
    dad_ht = 2*dad_a + dad_b

    # save both sibs haplotype indices
    sibs_ix[i] = mom_ht[0],dad_ht[0], mom_ht[1],dad_ht[1]


#reshape sibs_ix array so that it is one dimensional
#this contains the indices of the haplotype array in the order in which sibs haplotypes should be selected
sibs_ix = np.reshape(sibs_ix, (18000))

#reorder haplotypes in the previous generation in the order above
ht_sibs = ht[:,sibs_ix]

#save haplotype to file
np.savetxt(args.outpre + ".hap",
           ht_sibs,
           fmt="%u",
           delimiter=" ")

#write .legend file
chrom=vcf['variants/CHROM']
pos=vcf['variants/POS']
ref=np.full((nsnps), "A")
alt=np.full((nsnps), "T")

var_id = [str(a)+":"+str(b)+"_"+str(c)+"_"+str(d) for a,b,c,d in zip(chrom,pos,ref,alt)]

legend = np.column_stack((var_id,
                          vcf['variants/POS'],
                          vcf['variants/REF'],
                          vcf['variants/ALT']))

np.savetxt(args.outpre + ".legend",
           legend,
           header="id position a0 a1",
           delimiter=" ",
           fmt="%s",
          comments='')
