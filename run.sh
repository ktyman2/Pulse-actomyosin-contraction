#!/bin/bash
#SBATCH --job-name="pulse"
#SBATCH --output="run.%j.%N.out"
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --export=ALL
#SBATCH -t 48:00:00

ibrun -v ./main
