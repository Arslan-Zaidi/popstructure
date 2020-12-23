#!/bin/bash

vcf=${1}
out=${2}

bcftools view -i 'MAF>0.05' ${vcf} | gzip -> ${out}
