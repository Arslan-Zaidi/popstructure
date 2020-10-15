
for rep in {1..20}; do \

cat phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
> revisions/combined_sample/phenotypes/pheno_gridt100_ge_s9k.combined.${rep}.txt; \

tail -n+2 -q test/phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
>> revisions/combined_sample/phenotypes/pheno_gridt100_ge_s9k.combined.${rep}.txt; \

bash scripts/gwas/gwas.sh \
revisions/combined_sample/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.combined \
revisions/combined_sample/phenotypes/pheno_gridt100_ge_s9k.combined.${rep}.txt \
revisions/combined_sample/gwas_results/gwas_gridt100_combined.ge.${rep}.pcs0; \
done

gzip revisions/combined_sample/gwas_results/*.linear

for rep in {1..20}; do \
Rscript scripts/prs/clump_pcs0.R \
revisions/combined_sample/gwas_results/gwas_gridt100_combined.ge.${rep} \
smooth \
revisions/combined_sample/betas/est_effects.${rep}.smooth; \
done

for rep in {1..20}; do \
plink2 --pfile \
sibs/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.sib \
--score revisions/combined_sample/betas/est_effects.${rep}.smooth.nc.betas cols=scoresums \
--out revisions/combined_sample/prs/gridt100_prs_smooth.${rep}.nc \
--score-col-nums 3; \
done

for rep in {1..20}; do \
tail -n+2 -q revisions/combined_sample/prs/gridt100_prs_smooth.${rep}.nc.sscore | \
awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > revisions/combined_sample/prs/gridt100_prs_smooth.combined.all.nc.sscore
