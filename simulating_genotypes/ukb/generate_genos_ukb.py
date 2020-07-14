
#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--popfile","-p",dest="pop",help="file with info about latitude, longitude, sampling etc.",type=str,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
req_grp.add_argument("--migrate","-m",dest="migrate",help="migration rate",type=float,required=True)

parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
#parser.add_argument("--sampling",dest="scheme",help="sampling scheme (either uniform or weighted)",type=str,default="uniform",nargs="?")

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

print("loading data files")
#load adjacency matrix
#this will be the base for our migration matrix
adj_mat=np.genfromtxt("uk_nuts2_adj.txt",delimiter=",")

#load column names and rownames for matrix
mig_mat_id=pd.read_csv("uk_nuts2_adj_ids.txt",header=None,names=["deme_id"],sep="\t")

#multiply with scaling factor - the migration rate
mig_mat=args.migrate*adj_mat

#add columns and rows for ancient populations (with zeros)
mig_mat = np.vstack(( mig_mat, np.zeros( (3,35) )))
mig_mat = np.column_stack( (mig_mat, np.zeros((38,3)) ))

#load latitude and longitude for each deme
bplace_summary=pd.read_csv(args.pop,sep="\t")

#count number of individuals in each deme
bplace_summary=bplace_summary.groupby("deme").size().reset_index(name="ninds")

#merge migmat_ids with bplace_summary to preserve the same order as the migration matrix
bplace_summary2=mig_mat_id.merge(bplace_summary,left_on="deme_id",right_on="deme")

#number of demes
d=35

# #diploid sample size per deme
# ss_each=args.ss
#
# total_inds=d*args.ss
#
# if(args.scheme=="uniform"):
#     sample_sizes=[ss_each]*d
#
# if(args.scheme=="weighted"):
#     total_inds_ukb=np.sum(bplace_summary2['ninds'])
#     sample_sizes=(bplace_summary2['ninds']/total_inds_ukb) * total_inds
#     sample_sizes=[round(i) for i in sample_sizes]

sample_sizes = bplace_summary2['ninds']*2

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
    [msprime.MassMigration(time=100.1, source=i, destination=34, proportion=1.0) for i in range(6,12)],
    [msprime.MassMigration(time=100.1, source=i, destination=34, proportion=1.0) for i in range(12,18)],
    [msprime.MassMigration(time=100.1, source=i, destination=34, proportion=1.0) for i in range(18,24)],
    [msprime.MassMigration(time=100.1, source=i, destination=34, proportion=1.0) for i in range(24,30)],
    [msprime.MassMigration(time=100.1, source=i, destination=34, proportion=1.0) for i in range(30,34)],

    [msprime.MigrationRateChange(time=100.2,rate=0.004,matrix_index=(0,34))],
    [msprime.MigrationRateChange(time=100.2,rate=0.004,matrix_index=(34,0))],

    #ancient history
    #t 4500: Migrate lineages from 34 > 35 (WHG-south to Steppe) & from 0 > 35 (WHG-north to steppe)
    [msprime.MassMigration(time=4500, source=34, destination=35, proportion=0.2)],
    [msprime.MassMigration(time=4501, source=0, destination=35, proportion=0.5)],

    #t 7510: Migrate lineages from 34 > 38 (WHG-south > EF) & 0>37 (WHG-north > EF)
    [msprime.MassMigration(time=7510, source=34, destination=37, proportion=0.75)],
    [msprime.MassMigration(time=7511, source=0, destination=37, proportion=0.4)],

    #t 9000: Migrate lineages 35 > 36: steppe to EHG & WHG1 and WHG2 merge
    [msprime.MassMigration(time=9000, source=35, destination=36, proportion=0.5)],
    [msprime.MassMigration(time=9001, source=34, destination=0, proportion=1)],

    [msprime.MigrationRateChange(time = 9001, rate=0)],

    #t 25k: Migrate lineages 35 > 37, steppe & EF merge
    [msprime.MassMigration(time=25000, source=35, destination=37, proportion=1)],

    #t 30k: Migrate lineages 0 > 36: WHG and EHG merge
    [msprime.MassMigration(time=30000, source=0, destination=36, proportion=1)],

    #t 45k: Migrate lineages 36 > 37: HG and basal Eurasians merge
    [msprime.MassMigration(time=45000, source=36, destination=37, proportion=1)]]

demog = [item for sublist in demog_list for item in sublist]

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
