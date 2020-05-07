

#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
req_grp.add_argument("--migrate","-m",dest="migrate",help="migration rate",type=float,required=True)
args=parser.parse_args()

print(args)


import msprime
import numpy as np
import statistics
import math
import allel
import pandas as pd
import statsmodels.api as sm
from scipy import (stats,ndimage)

#number of demes in square
d=36

# dimension of square grid
dim=6

# define 2d grid to with deme identity
pmat=np.arange(0,d).reshape(dim,dim)

#define function to generate adjacency matrix
#arguments:
#m = migration rate in one direction
#nd = number of demes
def step_mig_mat(m,nd):
    #m is the uni-directional symmetric migration rate
    #NOTE: nd MUST be a squared number
    if(math.sqrt(nd).is_integer()==False):
        raise ValueError("nd must be a squared number (e.g. 4, 9, 16 ...) for the 2D model")
    else:
        nd2=int(math.sqrt(nd))
        #create matrix which will be used to determine which cells are adjacent in 2-dimensions
        #diagonals not considered for now but can be incorporated later if needed
        pmat=np.arange(0,nd).reshape(nd2,nd2)

        #create empty migration matrix to be filled in. This will be the output
        mmat=np.zeros(shape=[ (nd+3),(nd+3) ])

        #go through each cell in pmat and find out which cells are adjacent
        #first define functions to take care of corners and sides
        def contain(ix,max_ix):
            if ix<0:
                return(0)
            if ix>(max_ix-1):
                return(max_ix-1)
            else:
                return(ix)

        for ii in range(nd):
            center_ix=np.where(pmat==ii)
            top_ix=pmat[contain(center_ix[0]-1,nd2),contain(center_ix[1],nd2)]
            bottom_ix=pmat[contain(center_ix[0]+1,nd2),contain(center_ix[1],nd2)]
            left_ix=pmat[contain(center_ix[0],nd2),contain(center_ix[1]-1,nd2)]
            right_ix=pmat[contain(center_ix[0],nd2),contain(center_ix[1]+1,nd2)]

            mmat[ii,top_ix]=mmat[ii,bottom_ix]=mmat[ii,left_ix]=mmat[ii,right_ix]=m
            mmat[top_ix,ii]=mmat[bottom_ix,ii]=mmat[left_ix,ii]=mmat[right_ix,ii]=m

            mmat[ii,ii]=0

    return(mmat)

#generate migration matrix with migration rate provided by user
mig_mat=step_mig_mat(m=args.migrate,nd=d)

#set up the initial population
ss_each=args.ss*2
sample_sizes=[ss_each]*(d)

population_configurations = [
msprime.PopulationConfiguration(sample_size=k)
for k in sample_sizes]

population_configurations.extend([msprime.PopulationConfiguration(sample_size=0)]*3)

############ set up the demography

demog_list=[
    #change migration rate to 0 on the 100th generation
    [msprime.MigrationRateChange(time=100,rate=0)],
    #move lineages to the north (deme = 0) or south (deme = 35) to create N-S gradient
    [msprime.MassMigration(time=100, source=i, destination=0, proportion=1.0) for i in range(1,6)],
    [msprime.MassMigration(time=100, source=i, destination=0, proportion=0.8) for i in range(6,12)],
    [msprime.MassMigration(time=100, source=i, destination=0, proportion=0.6) for i in range(12,18)],
    [msprime.MassMigration(time=100, source=i, destination=0, proportion=0.4) for i in range(18,24)],
    [msprime.MassMigration(time=100, source=i, destination=0, proportion=0.2) for i in range(24,30)],
    [msprime.MassMigration(time=100.1, source=i, destination=35, proportion=1.0) for i in range(6,12)],
    [msprime.MassMigration(time=100.1, source=i, destination=35, proportion=1.0) for i in range(12,18)],
    [msprime.MassMigration(time=100.1, source=i, destination=35, proportion=1.0) for i in range(18,24)],
    [msprime.MassMigration(time=100.1, source=i, destination=35, proportion=1.0) for i in range(24,30)],
    [msprime.MassMigration(time=100.1, source=i, destination=35, proportion=1.0) for i in range(30,35)],

    #ancient history
    #t 4500: Migrate lineages from 35 > 36 (WHG-south to Steppe) & from 0 > 36 (WHG1 to steppe)
    [msprime.MassMigration(time=4500, source=35, destination=36, proportion=0.2)],
    [msprime.MassMigration(time=4501, source=0, destination=36, proportion=0.5)],

    #t 7510: Migrate lineages from 35 > 38 (WHG-south > EF) & 0>38 (WHG-n > EF)
    [msprime.MassMigration(time=7510, source=35, destination=38, proportion=0.75)],
    [msprime.MassMigration(time=7511, source=0, destination=38, proportion=0.4)],

    #t 9000: Migrate lineages 36 > 37: steppe to EF & WHG1 and WHG2 merge
    [msprime.MassMigration(time=9000, source=36, destination=37, proportion=0.5)],
    [msprime.MassMigration(time=9001, source=35, destination=0, proportion=1)],

    #t 25k: Migrate lineages 36 > 38, steppe & EF merge
    [msprime.MassMigration(time=25000, source=36, destination=38, proportion=1)],

    #t 30k: Migrate lineages 0 > 37: HG and basal Eurasians merge
    [msprime.MassMigration(time=30000, source=0, destination=37, proportion=1)],

    #t 45k: Migrate lineages 37 > 38: HG and basal Eurasians merge
    [msprime.MassMigration(time=45000, source=37, destination=38, proportion=1)]]

demog = [item for sublist in demog_list for item in sublist]

ts=msprime.simulate(Ne=10000,
                      population_configurations=population_configurations,
                      migration_matrix=mig_mat,
                      mutation_rate=args.mu,
                      recombination_rate=args.rho,
                      length=args.length,
                   demographic_events=demog)


print("writing genotype to vcf file")

with open(args.outpre+"_chr"+args.chr+".vcf","w") as vcf_file:
    ts.write_vcf(vcf_file,ploidy=2,contig_id=args.chr)

# #write pop file (deme ids for each individual)
#
# #write deme id/subpopulation for each individual
# deme_id=[[i]*ss_each for i in range(0,d)]
# #flatten
# deme_id=[item for sublist in deme_id for item in sublist]
#
# #write longitude and latitude for each individual
# bplace_x=[]
# bplace_y=[]
# for i in range(0,dim):
#     bplace_y.extend( [i] * ss_each * dim )
#     bplace_x.extend([item for item in range(0,dim) for i in range(ss_each)])
#
# #fid and iid
# fid=["msp_"+str(i) for i in range(0,(ss_each*d))]
# iid=["msp_"+str(i) for i in range(0,(ss_each*d))]
#
# popdf=pd.DataFrame({"FID":fid,
#                   "IID":iid,
#                   "POP":deme_id,
#                   "Latitude":bplace_y,
#                   "Longitude":bplace_x})
#
# popdf.to_csv(args.outpre+".pop",sep="\t",header=False,index=False)
