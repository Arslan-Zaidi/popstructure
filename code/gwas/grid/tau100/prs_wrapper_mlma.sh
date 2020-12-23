#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs_mlma.sh smooth cm ${rep}
bash scripts/prs/cal_prs_mlma.sh smooth re ${rep}

# bash scripts/prs/cal_prs_mlma.sh sharp cm ${rep}
# bash scripts/prs/cal_prs_mlma.sh sharp re ${rep}
