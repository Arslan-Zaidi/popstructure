
#script carries out matings both random and assortative (for deme) to generate sib pairs

import argparse
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--vcf", "-v", dest="vcf", help="vcf file. needed only for its header", type=str, required=True)
parser.add_argument("--iteration", "-i", dest="iteration", help="iteration", type=str, default=1, nargs="?")
args=parser.parse_args()

import pandas as pd
import numpy as np
import random
import allel

#read vcf header which contains list of sample names etc.
vcf_header = allel.read_vcf_headers(args.vcf)

#extract sample names from header
samples = vcf_header.samples

#get the deme id for each sample
demes = [i for i in range(0,36) for j in range(0,250)]



####### Make mate pairs
# (1)randomly across the entire grid

pairs=[] # this will store mom and dad's IDs

# list to store the population id for the sibling
#this will be randomly picked to be the pop of one of the parents
sibs_pop=[]

npairs = len(samples)//2

#list of samples to draw pairs from
samples_remaining = samples

i = 0
# loop over and make pairs until there are npairs pairs
while i < npairs:
    pairs.append(random.sample(samples,k=2))

    #add poplation information. Use one of the parents
    #this will later be used to generate environmental effects
    mom = pairs[i][0]
    mom_pop = [demes[j] for j in range(0,9000) if samples[j]==mom]
    sibs_pop.extend(mom_pop)

    samples_remaining=[item for item in samples_remaining if item not in pairs[i]]
    i=i+1

# generate IDs for siblings.
# adjacent even/odd numbers are assigned to sibling pairs
# this is not super important though as they are also identified by their parent information
sibs_id = ["i" + args.iteration + "1_"+ str(i) for i in range(0,9000,2)] + ["i"+ args.iteration + "_" + str(i) for i in range(0,9000,2)] #sibling IDs

# paste mum and dad's IDs together to create group ID
group = [i+","+j for i,j in pairs]
sex = [0]*9000 # sex - leave unknown

ped = np.column_stack( (sibs_id,
                        sibs_pop + sibs_pop,
                        group + group,
                        sex))

#save sample information to file - use this as .sample file
np.savetxt("test/sibs/genotypes/geno_complex_random_" + args.iteration + ".sample",ped,
           delimiter=" ",
          header="sample population group sex",
          fmt="%s",
          comments="")


####### Make mate pairs
# (2) assortatively, mates only find their partners within their own deme

pairs=[]   # this will store mom and dad's IDs
sibs_pop=[] # list where the population information will be stored for each sib

#loop over each deme (outer loop)
for i in range(0,36):
    # sample names from deme i
    samples_remaining = [samples[l] for l in range(0,9000) if demes[l]==i]
    npairs = len(samples_remaining)//2
    pairs_pop = []

    j = 0
    # loop over and make pairs until there are npairs pairs (inner loop)
    while j < npairs:

        #sibs population information - same for both parents so just one
        pairs_pop.append(random.sample(samples_remaining,k=2))
        mom_pop = [demes[l] for l in range(0,9000) if samples[l]==mom]
        sibs_pop.extend(mom_pop)

        samples_remaining=[item for item in samples_remaining if item not in pairs_pop[j]]

        j=j+1

    pairs.extend(pairs_pop)

#write to file
ped = np.column_stack( (sibs_id,
                        sibs_pop+sibs_pop,
                        group+group,
                        sex))

np.savetxt("test/sibs/genotypes/geno_complex_assort_" + args.iteration + ".sample",
            ped,
            delimiter=" ",
            header="sample population group sex",
            fmt="%s",
            comments="")
