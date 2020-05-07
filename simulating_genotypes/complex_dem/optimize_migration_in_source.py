
#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,required=True)
#req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--nreps","-n",dest="nreps",help="no. of reps",type=int,required=True)
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
#parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
#parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
req_grp.add_argument("--migrate","-m",dest="migrate",help="migration rate",type=float,required=True)
args=parser.parse_args()

print(args)

import msprime
import numpy as np
import math
import allel
import pandas as pd

# set up the initial population
ss_each=args.ss*2
sample_sizes=[ss_each]*5

population_configurations = [
msprime.PopulationConfiguration(sample_size=k)
for k in sample_sizes]


def fmig2(m, ndemes=5):
    migration_matrix = np.zeros( (ndemes,ndemes) )
    migration_matrix[0,1] = migration_matrix[1,0] = m
    return migration_matrix

#population_configurations.extend([msprime.PopulationConfiguration(sample_size=0)]*3)

migration_matrix = fmig2(args.migrate, ndemes=5)

demog_list=[

    #ancient history
    #t 4500: Migrate lineages from 35 > 36 (WHG-south to Steppe) & from 0 > 36 (WHG1 to steppe)
    [msprime.MassMigration(time=4500, source=1, destination=2, proportion=0.2)],
    [msprime.MassMigration(time=4501, source=0, destination=2, proportion=0.5)],

    #t 7510: Migrate lineages from 35 > 38 (WHG-south > EF) & 0>38 (WHG-n > EF)
    [msprime.MassMigration(time=7510, source=1, destination=4, proportion=0.75)],
    [msprime.MassMigration(time=7511, source=0, destination=4, proportion=0.4)],

    #t 9000: Migrate lineages 36 > 37: steppe to EF & WHG1 and WHG2 merge
    [msprime.MassMigration(time=9000, source=2, destination=3, proportion=0.5)],
    [msprime.MassMigration(time=9001, source=1, destination=0, proportion=1)],

    [msprime.MigrationRateChange(time = 9001, rate=0)],

    #t 25k: Migrate lineages 36 > 38, steppe & EF merge
    [msprime.MassMigration(time=25000, source=2, destination=4, proportion=1)],

    #t 30k: Migrate lineages 0 > 37: HG and basal Eurasians merge
    [msprime.MassMigration(time=30000, source=0, destination=3, proportion=1)],

    #t 45k: Migrate lineages 37 > 38: HG and basal Eurasians merge
    [msprime.MassMigration(time=45000, source=3, destination=4, proportion=1)]]

demog = [item for sublist in demog_list for item in sublist]

reps=msprime.simulate(Ne=10000,
                    population_configurations=population_configurations,
                    migration_matrix = migration_matrix,
                    mutation_rate=1e-08,
                    recombination_rate=1e-08,
                    length=args.length,
                    demographic_events=demog,
                    num_replicates = args.nreps)

fstmat = np.zeros((args.nreps,1))
subpops = np.array([list(range( ((i-1)*args.ss) , ((i-1)*args.ss)+args.ss)) for i in range(1,3)])

for i,ts in enumerate(reps):
    h = ts.genotype_matrix()
    daf=np.mean(h,axis=1)
    cm_ix=(daf > 0.05) & (daf < 0.95)

    h2 = allel.HaplotypeArray(h[cm_ix,:])
    g = h2.to_genotypes(ploidy=2)
    g2 = g.to_n_alt()
    g2_pruned_list = allel.locate_unlinked(g2, size=100, step=10, threshold=0.1)
    g3 = g[g2_pruned_list,:]


    a,b,c =allel.weir_cockerham_fst(g3, subpops)
    a = np.take(a,indices=1,axis=1)
    b = np.take(b,indices=1,axis=1)
    c = np.take(c,indices=1,axis=1)

    fst2 = np.sum(a)/( np.sum(a) + np.sum(b) + np.sum(c) )
    fstmat[i,0] = fst2

print(np.mean(fstmat))
