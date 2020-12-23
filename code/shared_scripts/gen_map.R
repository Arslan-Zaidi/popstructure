#!/usr/bin/env Rscript

args = commandArgs(TRUE)

if(length(args)!=2){
  stop("Usage: gen_map.R <legend file> <output file path>")
}

library(data.table)
library(dplyr)
library(rprojroot)

F = is_rstudio_project$make_fix_file()
legend.file.name = args[1]
output_file = args[2]

lfile = fread(F(legend.file.name))

lfile[, gdistance := position/1e6 ]
lfile = lfile[,.(position,gdistance)]
lfile[, cm_mb := 1]

lfile = lfile[,.(position, cm_mb, gdistance)]
colnames(lfile) = c("position","COMBINED_rate(cM/Mb)","Genetic_Map(cM)")

fwrite(lfile, F(output_file),
       col.names = TRUE,
       row.names = FALSE,
       quote = FALSE, sep="\t")
