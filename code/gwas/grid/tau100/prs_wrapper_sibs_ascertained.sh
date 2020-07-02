#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs_sibs_ascertained.sh smooth ${rep}

bash scripts/prs/cal_prs_sibs_ascertained.sh sharp ${rep}
