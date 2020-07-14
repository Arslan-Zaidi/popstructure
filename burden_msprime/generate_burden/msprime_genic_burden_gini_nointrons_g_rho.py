
# script generates variants in 100 genes of 39Kb each for 9k diploids
# gene burden is then calculated across all exons (refer to gene_ranges.txt file) for rare variants (MAF<0.001)
# then association tests are performed between the burden and one randomly sampled variant against
#the phenotype used is a sharp effect in a single deme (36 replicates for each deme separately)

#input:
# required:
#   1. recombination_rate
#   2. random_seed: random_seed for the simulation.
#      this will also be used later to look at the specific trees that generated potentially 'interesting'

# output files:
#   1. pvalues for burden tests for 100 genes for each of 36 demes (100x36 matrix)
#   2. dispersion for burden and variant (100x2)

import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

#required arguments
req_grp.add_argument("--seed","-x",dest="seed",help="random seed",type=int,required=True)
req_grp.add_argument("--length","-L",dest="length",help="length of chromosome (bp)",type=int,required=True)


parser.add_argument("--rho","-r",dest="rho",help="recombination rate (default=0.0)",type=float,default=0.0,nargs="?")
parser.add_argument("--Ne","-N",dest="npop",help="effective pop size (def:1e4)",type=int,default=10000,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--migrate","-m",dest="migrate",help="migration rate (def:0.05)",type=float,default=0.05,nargs="?")
parser.add_argument("--ndemes","-d",dest="ndemes",help="number of demes - must be a squared number (e.g. 25 or 100)",type=int,default=36,nargs="?")
parser.add_argument("--tmove","-t",dest="tmove",help="time (g) to panmixia. can be -9 (inf) or any positive integer",type=int,default=100,nargs="?")
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
    #tmove is the number of generations past which all lineages are moved into one deme.
    	#The is to reduce computational time when the no. of lineages << ndemes
        #also to mimic migration of an ancient population after which structure is established
        #set to 1000 generations by default

    sample_sizes=[ss_each]*d

    population_configurations = [
    msprime.PopulationConfiguration(sample_size=k)
    for k in sample_sizes]


    if tmove==-9:
         ts=msprime.simulate(Ne=args.npop,
                          population_configurations=population_configurations,
                          migration_matrix=mig_mat,
                          mutation_rate=args.mu,
                          recombination_rate=args.rho,
                          length=args.length,
                            num_replicates=100,
                            random_seed=args.seed)
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

#define vectors, which will be used to select the odd and even haplotypes of an individual
evens=range(0,ss*2*d,2)
odds=range(1,ss*2*d,2)

#get burden for each gene
burden=[]
variants=np.empty((100,1))
for j, tree_sequence in enumerate(ts):
    dosage=[]
    nvars=0
    for variant in tree_sequence.variants():
            daf=np.mean(variant.genotypes)
            if(daf<0.001):
                dosage.append(variant.genotypes[evens]+variant.genotypes[odds])
                nvars=nvars+1

    variants[j,0]=nvars
    #aggregate across all such variants and calculate burden for each individual
    #check if the dosage is nonzero first
    #found a situation where there were no rare variants in a gene (with rho=1e-08)
    #if dosage==0, append nan. we will deal with this later
    if(len(dosage)==0):
        burden.append(np.repeat(np.nan,ss*d))
    else:
        burden.append(np.sum(dosage,axis=0))

#load phenotype file
pheno=pd.read_csv("phenotypes/grid_d36_s250_e2.pheno",sep="\t")

#load genetic pcs
if(args.tmove==100):
    gpc_com=np.loadtxt("gpcs/genos_grid_d36_m0.05_s250_t100_gl200.common.ldpruned.eigenvec",
                           usecols=range(2,102))

    gpc_rare=np.loadtxt("gpcs/genos_grid_d36_m0.05_s250_t100_gl200.rare.ldpruned.eigenvec",
                            usecols=range(2,102))

if(args.tmove==-9):
    gpc_com=np.loadtxt("gpcs/genos_grid_d36_m0.07_s250_t-9_gl200.common.ldpruned.eigenvec",
                           usecols=range(2,102))

    gpc_rare=np.loadtxt("gpcs/genos_grid_d36_m0.07_s250_t-9_gl200.rare.ldpruned.eigenvec",
                            usecols=range(2,102))

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
ngenes=len(burden)
np_p1=np.empty((ngenes,8))
for i in range(0,ngenes):
        #without correction
        stat_result1=assoc(burden[i],pheno['smooth'],gpc_com,npcs=0)
        stat_result2=assoc(burden[i],pheno['sharp'],gpc_com,npcs=0)

        #correction using 100 common pcs
        stat_result3=assoc(burden[i],pheno['smooth'],gpc_com,npcs=100)
        stat_result4=assoc(burden[i],pheno['sharp'],gpc_com,npcs=100)

        #correction using 100 rare pcs
        stat_result5=assoc(burden[i],pheno['smooth'],gpc_rare,npcs=100)
        stat_result6=assoc(burden[i],pheno['sharp'],gpc_rare,npcs=100)

        np_p1[i,0]=args.seed
        np_p1[i,1]=i
        np_p1[i,2]=stat_result1[1]
        np_p1[i,3]=stat_result2[1]
        np_p1[i,4]=stat_result3[1]
        np_p1[i,5]=stat_result4[1]
        np_p1[i,6]=stat_result5[1]
        np_p1[i,7]=stat_result6[1]


#now calculating the gini curve for the burden in each gene
pd1=pd.DataFrame({'deme':pheno['deme'],
                 'bplace_x':pheno['longitude'],
                 'bplace_y':pheno['latitude']})


for i in range(0,100):
    key1='burden_'+str(i)
    pd1[key1]=burden[i]

print("calculating gini curve")

##Calculate no. of demes the variant appears in
gini_b=np.empty((100,36))
bsum=[]
for k in range(0,100):
    key1='burden_'+str(k)

    x=pd1.groupby(['deme','bplace_x','bplace_y'],as_index=False)[key1].sum()
    x=x[key1]
    bsum.append(np.sum(x))

    x2=np.cumsum(np.sort(x))/np.sum(x)
    gini_b[k,:]=x2


final_output=np.column_stack((np.repeat(args.seed,100),np.arange(100),variants,bsum,gini_b))


print("final output, hang tight")
#output all files

#output pvalue for association test
np.savetxt("output/lvburden/association/t"+str(args.tmove)+"/pburden_x"+str(args.seed)+"_l"+str(args.length)+".txt",np_p1,delimiter=",",fmt='%6.6g')

#output dispersion
np.savetxt("output/lvburden/gini/t"+str(args.tmove)+"/pgini_x"+str(args.seed)+"_l"+str(args.length)+".txt",
           final_output,
           delimiter=",",
           fmt='%3.3g')
