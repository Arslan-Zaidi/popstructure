

#script generates simulated genotype data
#the primary motviation to write this is to generate genotype data (a time-consuming step), which can then be used for downstream analyses
import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="total sample size",type=int,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
parser.add_argument("--Ne","-N",dest="npop",help="effective pop size (def:1e4)",type=int,default=10000,nargs="?")
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")

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

print("simulating genotypes under the coalescent")

ts = msprime.simulate(Ne = args.npop,
                    mutation_rate = args.mu,
                    recombination_rate = args.rho,
                    length = args.length,
                    sample_size = args.ss)

print("writing genotype to vcf file")

with open(args.outpre+"_chr"+args.chr+".vcf","w") as vcf_file:
    ts.write_vcf(vcf_file,ploidy=2,contig_id=args.chr)
