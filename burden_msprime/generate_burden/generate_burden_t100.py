
# script generates variants in 100 genes of 50Kb each for 9k diploids
# gene burden is then calculated across all exons (refer to gene_ranges.txt file) for rare variants (MAF<0.001)
# then association tests are performed between the burden and one randomly sampled variant against
#the phenotype used is a sharp effect in a single deme (36 replicates for each deme separately)

#reference for gene structure:
#https://link.springer.com/article/10.1186/s13104-019-4343-8/tables/2

import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

#required arguments
req_grp.add_argument("--rho","-r",dest="rho",help="recombination rate",type=float,required=True)
req_grp.add_argument("--seed","-x",dest="seed",help="random seed",type=int,required=True)
req_grp.add_argument("--migrate","-m",dest="migrate",help="migration rate (def:0.05)",type=float,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output_prefix",type=str,required=True)

#preset arguments
parser.add_argument("--Ne","-N",dest="npop",help="effective pop size (def:1e4)",type=int,default=10000,nargs="?")
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:50kb)",type=int,default=50000,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--ndemes","-d",dest="ndemes",help="number of demes - must be a squared number (e.g. 25 or 100)",type=int,default=36,nargs="?")
parser.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,default=250,nargs="?")

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


#number of demes
d=args.ndemes

# dimension of square grid
dim=int(np.sqrt(d))

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
mig_mat=step_mig_mat(m=0.05,nd=d)

#diploid sample size within each deme
ss=args.ss

#number of haplotypes
nhaps=ss*2

##### define function to simulate genotypes under a stepping stone migration model
def step_geno(ss_each=ss*2,tmove=100):
    #N is the population size for each deme
    #ss_each is the haploid sample size for each deme
    #l is the length of the chromosome

    sample_sizes=[ss_each]*d

    population_configurations = [
    msprime.PopulationConfiguration(sample_size=k)
    for k in sample_sizes]

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


    ts=msprime.simulate(Ne=args.npop,
                          population_configurations=population_configurations,
                          migration_matrix=mig_mat,
                          mutation_rate=args.mu,
                          recombination_rate=args.rho,
                          length=args.length,
                       demographic_events=demog,
                       num_replicates=100,
                       random_seed=args.seed)

    return(ts)

print("simulating genealogies")
#simulate!
ts=step_geno(ss_each=ss*2,tmove=args.tmove)

print("calculating burden for each gene")
#if recombination rate == 0, calculate across entire gene (simulate shorter gene)
if args.rho == 0:
    #define vectors, which will be used to select the odd and even haplotypes of an individual
    evens=range(0,ss*2*d,2)
    odds=range(1,ss*2*d,2)

    #get burden and sample a single rare variant from each gene
    burden = np.empty((100,ss*args.ndemes))
    nvariants = np.zeros((100,1))
    for j, tree_sequence in enumerate(ts):
        dosage=[]
        variants=[]
        for variant in tree_sequence.variants():
                daf=np.mean(variant.genotypes)
                if(daf<0.001):
                    dosage.append(variant.genotypes[evens]+variant.genotypes[odds])
                    variants.append(1)

        nvariants[j,:]=np.sum(variants)
        #aggregate across all such variants and calculate burden for each individual
        #check if the dosage is nonzero first
        #found a situation where there were no rare variants in a gene (with rho=1e-08)
        #if dosage==0, append nan. we will deal with this later
        if(len(dosage)==0):
            burden[j,:] = np.repeat(np.nan,ss*d)
        else:
            burden[j,:] = np.sum(dosage,axis=0)

if args.rho>0:
    #describe the exon structure of the genes
    #we don't need the introns or the UTRs as we did in SLIM
    gene_ranges=[(0,160),(7098,7258),(14196,14356),(21294,21454),(28392,28552),(35490,35650),(42588,42748),(49686,49846)]

    #define vectors, which will be used to select the odd and even haplotypes of an individual
    evens=range(0,ss*2*d,2)
    odds=range(1,ss*2*d,2)

    #get burden and sample a single rare variant from each gene
    burden = np.empty((100,ss*args.ndemes))
    nvariants = np.zeros((100,1))
    for j, tree_sequence in enumerate(ts):
        dosage=[]
        variants=[]
        for variant in tree_sequence.variants():
            if any(lower<=variant.site.position<=upper for (lower,upper) in gene_ranges):
                daf=np.mean(variant.genotypes)
                if(daf<0.001):
                    dosage.append(variant.genotypes[evens]+variant.genotypes[odds])
                    variants.append(1)

        nvariants[j,:]=np.sum(variants)
        #aggregate across all such variants and calculate burden for each individual
        #check if the dosage is nonzero first
        #found a situation where there were no rare variants in a gene (with rho=1e-08)
        #if dosage==0, append nan. we will deal with this later
        if(len(dosage)==0):
            burden[j,:] = np.repeat(np.nan,ss*d)
        else:
            burden[j,:] = np.sum(dosage,axis=0)

burden=np.column_stack((np.repeat(args.seed,100),
                        np.arange(0,100),
                        burden))

#save burden to compressed file for future use
np.savez_compressed(args.outpre + '.haps.npz', burden)

#write number of variants to file
np.savetxt(args.outpre + '_nvariants.txt',nvariants,delimiter=",",fmt='%i')
