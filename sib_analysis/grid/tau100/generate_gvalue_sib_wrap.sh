
#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

effects_i=${1}

for j in {1..20}; do \
python generate_gvalue_sib.py --chr ${j} --effect_i ${effects_i}; \
echo $j; \
done

# for j in {1..20}; do \
# tail -n+2 -q gvalues/genos_gridt100_l1e7_ss750_m0.05_chr${j}_e${effects_i}_sibs_assort.sscore | awk -v j=$j -v OFS="\t" '{print j,$0}'; \
# done
