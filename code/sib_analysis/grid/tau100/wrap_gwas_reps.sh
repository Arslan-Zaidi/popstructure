
#!/bin/bash

effects_rep=${1}

mkdir -p gwas_results/e${effects_rep}

for i in {1..20}; \
do bsub -o gwas_logs/gwas_i${i}.e${effects_rep}.log \
-e gwas_logs/gwas_i${i}.e${effects_rep}.err \
-M 30000 -R "rusage[mem=30000]" \
bash gwas_sib_grid_wrapper.sh ${i} ${effects_rep} 9000; done
