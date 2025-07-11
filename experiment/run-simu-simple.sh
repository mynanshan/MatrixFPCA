#!/bin/bash
#SBATCH --job-name=simu-simple-matfpca
#SBATCH --account=def-jiguocao
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1:00:00
#SBATCH --output="output_files/simu-simplu-%j.out"

module load r

n=$1:-10
p=$2:-10
rt=$3:-3

Rscript experiment/simu-simple-full.R -n $n \
                               -p $p \
                               -q $p \
                               --rank_time $rt \
                               --compare_methods 0