
#script carries out matings both random and assortative (for deme) to generate sib pairs

import argparse
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--vcf", "-v", dest="vcf", help="vcf file. needed only for its header", type=str, required=True)
parser.add_argument("--iteration", "-i", dest="iteration", help="iteration", type=str, default=1, nargs="?")
req_grp.add_argument("--outpre", "-o", dest="outpre", help="output prefix", type=str, required=True)
args=parser.parse_args()

import numpy as np
import random
import allel

random.seed(args.iteration)

#read vcf header which contains list of sample names etc.
vcf_header = allel.read_vcf_headers(args.vcf)

#extract sample names from header
samples = vcf_header.samples

#get the deme id for each sample
demes = [i for i in range(0,36) for j in range(0,250)]



# ####### Make mate pairs
# # (1)randomly across the entire grid
#
# pairs=[] # this will store mom and dad's IDs
#
# # list to store the population id for the sibling
# #this will be randomly picked to be the pop of one of the parents
# sibs_pop=[]
#
# npairs = len(samples)//2
#
# #list of samples to draw pairs from
# samples_remaining = samples
#
# i = 0
# # loop over and make pairs until there are npairs pairs
# while i < npairs:
#     pairs.append(random.sample(samples,k=2))
#
#     #add poplation information. Use one of the parents
#     #this will later be used to generate environmental effects
#     mom = pairs[i][0]
#     mom_pop = demes[samples.index(mom)]
#     sibs_pop.append(mom_pop)
#
#     samples_remaining=[item for item in samples_remaining if item not in pairs[i]]
#     i=i+1
#
# #generate ID for each sibling pair
# pair_id = ["f" + str(args.iteration) + "_"+ str(i) for i in range(0,4500)] #sibling pair ID
#
# # paste mum and dad's IDs together to create group ID
# group = [i+","+j for i,j in pairs]
# sex = [0]*4500 # sex - leave unknown
#
# #get longitude and latitude for each deme
# lon = [ i%6 for i in sibs_pop]
# lat = [int(np.floor(i/6)) for i in sibs_pop]
#
# # I'm going to generate sibling genotypes in alternating columns
# # odd no.s for sibling 1 and even numbers for sibling 2
# sib1 = np.column_stack( (pair_id,
#                          ["i" + str(args.iteration) + "_" + str(i) for i in range(0,9000,2)],
#                         sibs_pop,
#                         lon,
#                         lat,
#                         group) )
#
# sib2 = np.column_stack( (pair_id,
#                          ["i" + str(args.iteration) + "_" + str(i) for i in range(1,9000,2)],
#                         sibs_pop,
#                         lon,
#                         lat,
#                         group) )
#
# sib_ped = np.empty((sib1.shape[0]*2,sib1.shape[1]), dtype=sib1.dtype)
#
# sib_ped[0::2] = sib1
# sib_ped[1::2] = sib2
#
#
# np.savetxt(args.outpre + "_random_" + str(args.iteration) + ".sample",
#           sib_ped,
#           delimiter=" ",
#           header="fid iid population longitude latitude group",
#           fmt="%s",
#           comments="")


####### Make mate pairs
# (2) assortatively, mates only find their partners within their own deme

pairs=[]   # this will store mom and dad's IDs
sibs_pop=[] # list where the population information will be stored for each sib

npairs = len(samples)//2
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
        mom = pairs_pop[j][0]
        mom_pop = demes[samples.index(mom)]
        sibs_pop.append(mom_pop)

        samples_remaining=[item for item in samples_remaining if item not in pairs_pop[j]]

        j=j+1

    pairs.extend(pairs_pop)

# #generate ID for each sibling pair
pair_id = ["f" + str(args.iteration) + "_"+ str(i) for i in range(0,4500)] #sibling pair ID

#new parents
group = [i+","+j for i,j in pairs]

#get longitude and latitude for each deme
lon = [ i%6 for i in sibs_pop]
lat = [int(np.floor(i/6)) for i in sibs_pop]

sib1 = np.column_stack( (pair_id,
                         ["i" + str(args.iteration) + "_" + str(i) for i in range(0,9000,2)],
                        sibs_pop,
                        lon,
                        lat,
                        group) )

sib2 = np.column_stack( (pair_id,
                         ["i" + str(args.iteration) + "_" + str(i) for i in range(1,9000,2)],
                        sibs_pop,
                        lon,
                        lat,
                        group) )

sib_ped = np.empty((sib1.shape[0]*2,sib1.shape[1]), dtype=sib1.dtype)

sib_ped[0::2] = sib1
sib_ped[1::2] = sib2

#write to file
np.savetxt(args.outpre + "_assort_" + args.iteration + ".sample",
            sib_ped,
            delimiter=" ",
            header="fid iid population longitude latitude group",
            fmt="%s",
            comments="")
