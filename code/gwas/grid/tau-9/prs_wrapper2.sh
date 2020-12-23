#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs2.sh smooth ${rep}

#bash scripts/prs/cal_prs2.sh smooth_long ${rep}

bash scripts/prs/cal_prs2.sh sharp ${rep}

#bash scripts/prs/cal_prs2.sh grandom ${rep}
