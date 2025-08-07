#!/bin/bash
#SBATCH --job-name=wofost_runs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --array=0-26
#SBATCH --output=wofost_job_%A_%a.out
#SBATCH --error=wofost_job_%A_%a.err

# Load required software modules
module load legacy
module load netcdf

# Define tsum and sow arrays
tsum_pairs=("e1e1" "e1a1" "e1l1" "a1e1" "a1a1" "a1l1" "l1e1" "l1a1" "l1l1")
sow_vars=("e1" "a1" "l1")

# Calculate indices
index=$SLURM_ARRAY_TASK_ID
tsum_idx=$((index / 3))
sow_idx=$((index % 3))

# Select tsum and sow values
tsum=${tsum_pairs[$tsum_idx]}
sow=${sow_vars[$sow_idx]}

# Execute the WOFOST command
/usr/bin/time -v -- ./wofost list.txt meteolist.txt all_griddata.nc "avg_tsum1_${tsum}" "avg_tsum2_${tsum}" "sow_${sow}" --use-potential-nutrients
