#!/bin/bash
#PBS -l nodes=5:ppn=1
#PBS -l walltime=24:00:00
#PBS -l mem=8gb

ml purge; ml R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

cd $PBS_O_WORKDIR
# Rscript run.r


