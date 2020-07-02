#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs_sibs.sh smooth random ${rep}
bash scripts/prs/cal_prs_sibs.sh smooth assort ${rep}

#bash scripts/prs/cal_prs_sibs.sh smooth_long ${rep}

bash scripts/prs/cal_prs_sibs.sh sharp random ${rep}
bash scripts/prs/cal_prs_sibs.sh sharp assort ${rep}


#bash scripts/prs/cal_prs_sibs.sh grandom ${rep}
