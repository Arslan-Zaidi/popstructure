

#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
parser.add_argument("--Ne","-N",dest="npop",help="effective pop size (def:1e4)",type=int,default=10000,nargs="?")
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--migrate","-m",dest="migrate",help="migration rate (def:0.05)",type=float,default=0.05,nargs="?")
parser.add_argument("--ndemes","-d",dest="ndemes",help="number of demes - must be a squared number (e.g. 25 or 100",type=int,default=36,nargs="?")
parser.add_argument("--tmove","-t",dest="tmove",help="time (g) to panmixia. can be -9 (inf) or any positive integer",type=int,default=100,nargs="?")
args=parser.parse_args()

print(args)

import msprime
import numpy as np
import statistics
import math
import allel
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd
from scipy import (stats,ndimage)

#number of demes
d=args.ndemes

#dimension of square grid
dim=int(np.sqrt(d))

#define 2d grid to with deme identity
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
        mmat=np.zeros(shape=[nd,nd])

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

#diploid sample size within each deme
ss=args.ss

##### define function to simulate genotypes under a stepping stone migration model
def step_geno(N=1e4,l=1e7,ss_each=2*ss,tmove=1000,mmat=mig_mat):
    #N is the population size for each deme
    #ss_each is the haploid sample size for each deme
    #l is the length of the chromosome
    #tmove is the number of generations past which all lineages are moved into one deme.
    	#The is to reduce computational time when the no. of lineages << ndemes
        #also to mimic migration of an ancient population after which structure is established
        #set to 1000 generations by default

    sample_sizes=[ss_each]*d

    population_configurations = [
    msprime.PopulationConfiguration(sample_size=k)
    for k in sample_sizes]


    if tmove==-9:
         ts=msprime.simulate(Ne=N,
                          population_configurations=population_configurations,
                          migration_matrix=mmat,
                          mutation_rate=args.mu,
                          recombination_rate=args.rho,
                          length=l)
    else:
        #specify demographic event - move all lineages to one population after tmove generations
        demog=[
            msprime.MassMigration(
                time=tmove,
                source=i,
                destination=d-1,
                proportion=1.0) for i in range(d-1)]

        demog.append(#change migration rate among demes to be 0
            msprime.MigrationRateChange(
                time=tmove,
                rate=0))


        ts=msprime.simulate(Ne=N,
                              population_configurations=population_configurations,
                              migration_matrix=mmat,
                              mutation_rate=args.mu,
                              recombination_rate=args.rho,
                              length=l,
                           demographic_events=demog)

    return(ts)

print("simulating genotypes under demographic model")

#simulate!
ts=step_geno(N=args.npop,ss_each=2*ss,l=args.length,tmove=args.tmove)

print("writing genotype to vcf file")

with open(args.outpre+"_chr"+args.chr+".vcf","w") as vcf_file:
    ts.write_vcf(vcf_file,ploidy=2,contig_id=args.chr)

#write pop file (deme ids for each individual)

#write deme id/subpopulation for each individual
deme_id=[[i]*ss for i in range(0,d)]
#flatten
deme_id=[item for sublist in deme_id for item in sublist]

#write longitude and latitude for each individual
bplace_x=[]
bplace_y=[]
for i in range(0,dim):
    bplace_y.extend( [i] * ss * dim )
    bplace_x.extend([item for item in range(0,dim) for i in range(ss)])

#fid and iid
fid=["msp_"+str(i) for i in range(0,(ss*d))]
iid=["msp_"+str(i) for i in range(0,(ss*d))]

popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id,
                  "Latitude":bplace_y,
                  "Longitude":bplace_x})

popdf.to_csv(args.outpre+".pop",sep="\t",header=False,index=False)
