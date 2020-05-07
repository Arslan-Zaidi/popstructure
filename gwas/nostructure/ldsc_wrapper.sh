
conda activate rpkgs

phenotype=${1}

Rscript ldsc_regression.R \
train/ldsc/chr1_20.ldscore.gz \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.2.pcs0.${phenotype}.glm.linear.gz \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.frq.afreq \
train/ldsc/ldsc_nostr_pcs0_${phenotype}

Rscript ldsc_regression.R \
train/ldsc/chr1_20.ldscore.gz \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.2.cm.${phenotype}.glm.linear.gz \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.frq.afreq \
train/ldsc/ldsc_nostr_cm_${phenotype}

Rscript ldsc_regression.R \
train/ldsc/chr1_20.ldscore.gz \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.2.re.${phenotype}.glm.linear.gz \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.frq.afreq \
train/ldsc/ldsc_nostr_re_${phenotype}

Rscript ldsc_regression.R \
train/ldsc/chr1_20.ldscore.gz \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.2.cmre.${phenotype}.glm.linear.gz \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.frq.afreq \
train/ldsc/ldsc_nostr_cmre_${phenotype}
