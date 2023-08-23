#!/bin/bash

#SBATCH --partition=fasse
#SBATCH --qos=normal
#SBATCH --account=dominici_lab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --ntasks-per-node=1
#SBATCH --job-name GPCERF_sim
#SBATCH --output GPCERF_sim.out
#SBATCH --mem=184GB
#SBATCH --time=02-00:00
#SBATCH --mail-user=mcork@g.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/9.3.0-fasrc01 R/4.0.5-fasrc02

R CMD BATCH --no-save "--args $1 $2" small_simulation.R Logs/data_app_fasse_$2.Rout
