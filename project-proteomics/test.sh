#!/bin/sh

#SBATCH --time 00:10:00
#SBATCH -o hello_bioc6242_%j.out
#SBATCH -e hello_bioc6242_%j.err
#SBATCH -p defq

echo "Hello World"

module load R
Rscript /groups/bioc6242/test.R
