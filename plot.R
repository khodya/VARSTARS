#!/usr/bin/Rscript

# Plotting script for HD52721 lightcurve modeling

# Load observation data
obs <-  read.table("obs/pmm_all_1610155_15")

# Load model data
system("./a.out")
calc <- read.table("output.txt")



