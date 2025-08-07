#!/bin/bash

#SBATCH --job-name=wofost_run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=wofost_job_%j.out
#SBATCH --error=wofost_job_%j.err

# Load required software modules
module load legacy
module load netcdf

# Execute the WOFOST command
./wofost list.txt meteolist.txt all_griddata.nc avg_tsum1_e1e1 avg_tsum2_e1e1 sow_e1
