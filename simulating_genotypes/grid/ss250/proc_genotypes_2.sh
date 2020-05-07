#!/bin/bash

### script converts vcf files for separate chromosomes to plink format and then merges the files

#exit and print usage if arguments not provided
[ $# -ne 4 ] && { echo "Usage: $0 <tau {100,-9}> <chr_number> <sample_size> <migration_rate>"; exit 1; }

tau=${1} # time to panmixia
nchr=${2} # no. of chromosomes to simulate
ss=${3} # sample size for EACH deme. NOT the total sample size
m=${4}

###### merging separate chromosomes into one file ####
cd genotypes/tau${tau}

echo "creating rslist in bim file"
#create snp id for each variant
#define snp name as combination of chromosome number and variant number
#e.g. rs2_300 refers to the 300th variant on the 2nd chromosome
for i in $(seq 1 $nchr); do \
awk '{print $1,"rs"$1"_"$4, $3,$4,$5,$6}' genos_grid_d36_m${m}_s${ss}_t${tau}_chr${i}.bim > tmp && mv tmp genos_grid_d36_m${m}_s${ss}_t${tau}_chr${i}.bim;
done

echo "creating list of files to be merged"
#making list of files to be merged
for i in $(seq 2 $nchr); do \
echo genos_grid_d36_m${m}_s${ss}_t${tau}_chr${i}.bed genos_grid_d36_m${m}_s${ss}_t${tau}_chr${i}.bim genos_grid_d36_m${m}_s${ss}_t${tau}_chr${i}.fam;
done > mergelist.txt

echo "merging files"
#calculate total length of genome (in Mb) _gl
chrlength=`expr $nchr \* 10`
#merge all chromosomes into one file
plink --bfile genos_grid_d36_m${m}_s${ss}_t${tau}_chr1 \
--merge-list mergelist.txt \
--make-bed \
--out genos_grid_d36_m${m}_s${ss}_t${tau}_gl${chrlength}

mkdir -p phenotypes

echo "cleaning up"
rm *file
rm *.nosex
#rm genos_grid_d36_m1e-03_s${ss}_t${tau}_chr*
rm mergelist.txt
# for i in $(seq 1 $nchr); do \
# #rm genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${i}.bed genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${i}.bim genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${i}.fam;
# done
