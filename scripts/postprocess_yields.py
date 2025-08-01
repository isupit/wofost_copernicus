import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Define the scenarios (matching the SLURM script)
tsum_pairs = ["e1e1", "e1a1", "e1l1", "a1e1", "a1a1", "a1l1", "l1e1", "l1a1", "l1l1"]
sow_vars = ["e1", "a1", "l1"]  # Corresponding to sow_e1, sow_a1, sow_l1

# Collect mean yield averages
data = []
for tsum in tsum_pairs:
    for sow in sow_vars:
        filename = f"wofost_results_{tsum}_sow_{sow}.nc"
        try:
            ds = nc.Dataset(filename, 'r')
            yield_avg = ds.variables['Yield_Average'][:]
            yield_avg_mean = np.nanmean(yield_avg)  # Global mean, ignoring NaN
            data.append({'tsum_pair': tsum, 'sowing_date': f"sow_{sow}", 'mean_yield': yield_avg_mean})
            ds.close()
        except FileNotFoundError:
            print(f"Warning: File {filename} not found, skipping.")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

# Create DataFrame for easy plotting
df = pd.DataFrame(data)
if df.empty:
    print("No data found. Ensure NetCDF files are present.")
    exit(1)

# Pivot for grouped bar chart (tsum_pair as index, sowing_date as columns)
df_pivot = df.pivot(index='tsum_pair', columns='sowing_date', values='mean_yield')

# Plot grouped bar chart
fig, ax = plt.subplots(figsize=(12, 6))
bar_width = 0.25
x = np.arange(len(tsum_pairs))

# Plot bars for each sowing date
ax.bar(x - bar_width, df_pivot['sow_e1'], width=bar_width, label='sow_e1 (early)', color='skyblue')
ax.bar(x, df_pivot['sow_a1'], width=bar_width, label='sow_a1 (average)', color='lightgreen')
ax.bar(x + bar_width, df_pivot['sow_l1'], width=bar_width, label='sow_l1 (late)', color='salmon')

# Labels and styling
ax.set_xlabel('Temperature Sum Pairs')
ax.set_ylabel('Mean Yield Average (kg/ha)')
ax.set_title('Mean Yield Averages Across Scenarios')
ax.set_xticks(x)
ax.set_xticklabels(tsum_pairs, rotation=45, ha='right')
ax.legend(title='Sowing Date')
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig('yield_averages_barchart.png')
print("Bar chart saved as 'yield_averages_barchart.png'")
