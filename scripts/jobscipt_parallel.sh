#!/bin/bash
#SBATCH --job-name=wofost_runs
#SBATCH --ntasks=1
#SBATCH --array=0-26  # Limit to 4 concurrent tasks to avoid OOM
#SBATCH --mem=4G  # Request 4 GB per job
#SBATCH --output=slurm-%A_%a.out
tsum_pairs=("e1e1" "e1a1" "e1l1" "a1e1" "a1a1" "a1l1" "l1e1" "l1a1" "l1l1")
sow_vars=("e1" "a1" "l1")
index=$SLURM_ARRAY_TASK_ID
tsum_idx=$((index / 3))
sow_idx=$((index % 3))
tsum=${tsum_pairs[$tsum_idx]}
sow=${sow_vars[$sow_idx]}
# Log memory usage for tuning
/usr/bin/time -v ./wofost list.txt meteolist.txt all_griddata_cropped.nc "avg_tsum1_${tsum}" "avg_tsum2_${tsum}" "sow_${sow}"
