
#this script carries out sibling GWAS using phenotypes and haplotypes from siblings

import argparse
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--haplotypes", "-g", dest="ht", help="haplotype array in npz format", type=str, required=True)
req_grp.add_argument("--phenotypes", "-p", dest="phenotypes", help="file with phenotypes per person", type=str, required=True)
req_grp.add_argument("--outpre", "-o", dest="outpre", help="prefix for output file - should be informative about mating", type=str, required=True)
req_grp.add_argument("--legend", "-l", dest="legend", help="file with snp information", type=str, required=True)
#parser.add_argument("--iteration", "-i", dest="iteration", help="iteration", type=str, default="1", nargs="?")
args=parser.parse_args()

import numpy as np
import statistics
import math
import allel
import statsmodels.api as sm
import pandas as pd
from scipy import stats

#read compressed haplotype array
ht = np.load(args.ht)
ht = ht['arr_0']

nsnps = ht.shape[0]
ninds = ht.shape[1]//2

#read legend file - snps info
snps = pd.read_csv(args.legend,
                   delim_whitespace=True)

#calculate allele frequency for each variant and filter out
daf = np.mean(ht, axis = 1)

#keep SNPs that are polymorphic
polymorphic_ix = np.where((daf>0) & (daf<1))[0]

ht = ht[polymorphic_ix,:]
snps = snps.iloc[polymorphic_ix]

#convert to scikit.allel genotype format
ht = allel.HaplotypeArray(ht)
gt = ht.to_genotypes(ploidy=2)

#convert gt to alternate allele count
gt = gt.to_n_alt()

#calculate the difference in genotype value between snps
#odd_sibs = gt[:,0::2] # odd sibs
#even_sibs = gt[:,1::2] # even sibs
gtd = gt[:,0::2] - gt[:,1::2]

#define function to test association between phenotype and a single variant
#define function to test association between phenotype and a single variant
def assoc(gmat,response):
    #gmat is the genotype matrix
    #response is the response variable
    #variant index is the SNP index
    #predictors=gmat[variant_index,:]
    if(np.mean(gmat)!=0):
        predictors=sm.add_constant(gmat)
        model=sm.OLS(response,predictors).fit()
        beta=model.params[1]
        se=model.bse[1]
        pval=model.pvalues[1]
        #calculate chisq (for calculation of lambda)
        chi=(beta/se)**2
        chi_p=1-stats.chi2.cdf(chi,1)
        results=[beta,se,pval,chi,chi_p]
    else:
        results=[np.nan]*5

    return results

#load phenotype file
phenotypes = pd.read_csv(args.phenotypes,
                        sep="\t")

# give a number to each sibling within each pair
phenotypes['sibling']=phenotypes.groupby('FID').cumcount().add(1)

# the phenotype file needs to be sorted to make sure it aligns with the genotype file
# to ensure this, add iteration and family no.
arr1 = np.array([i.split("_") for i in phenotypes.FID])
arr1[:,0]=np.array([i[1] for i in arr1[:,0]])
arr1=arr1.astype(int)

phenotypes['iteration_ix'] = arr1[:,0]
phenotypes['family_ix'] = arr1[:,1]

phenotypes=phenotypes.sort_values(by=['iteration_ix','family_ix'])

#pivot phenotype table so the phenotypes of siblings are side-by-side
#calculate difference in phenotype between sibs
phenotypes_diff=np.empty( ((ninds//2) , 4) )

def gwas(pheno_name):
    pivoted_df = pd.pivot_table(phenotypes,
                                index=['iteration_ix','family_ix'],
                                columns="sibling",
                                values=pheno_name)
    phenotypes_diff = pivoted_df.iloc[:,0] - pivoted_df.iloc[:,1]
    gwas_results = np.apply_along_axis(assoc,1,gtd,phenotypes_diff)
    gwas_results = np.column_stack((snps,gwas_results))
    gwas_results = pd.DataFrame(gwas_results,
                                columns=["SNP","POS","REF","ALT","BETA","SE","P","CHI","CHI_P"])
    return gwas_results


#grandom = gwas("grandom")
smooth  = gwas("smooth")
#smooth_long = gwas("smooth_long")
sharp = gwas("sharp")

# grandom.to_csv(args.outpre+".grandom.glm.linear",
#                 header=True,
#                 index=False,
#                 sep="\t")

smooth.to_csv(args.outpre+".smooth.glm.linear",
                header=True,
                index=False,
                sep="\t")

# smooth_long.to_csv(args.outpre+".smooth_long.glm.linear",
#                 header=True,
#                 index=False,
#                 sep="\t")
#
sharp.to_csv(args.outpre+".sharp.glm.linear",
                header=True,
                index=False,
                sep="\t")


#write to file
# np.savetxt(args.outpre+ "glm.linear",
#             gwas_results,
#             delimiter=" ",
#             header="SNP position ref alt  group",
#             fmt="%s",
#             comments="")
