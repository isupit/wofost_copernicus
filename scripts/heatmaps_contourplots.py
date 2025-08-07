import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import seaborn as sns
from matplotlib import patheffects as pe

# --- Configuration ---
# Define the directories containing your WOFOST output NetCDF files
data_dirs = {
    "original": "/lustre/scratch/WUR/ESG/datad002/iwan_utam/wofost_copernicus/output",
    "evtra": "/lustre/scratch/WUR/ESG/datad002/iwan_utam/wofost_copernicus/output_potential_evtra",
    "nutrients": "/lustre/scratch/WUR/ESG/datad002/iwan_utam/wofost_copernicus/output_potential_nutrients",
    "evtra_nutrients": "/lustre/scratch/WUR/ESG/datad002/iwan_utam/wofost_copernicus/output_potential_evtra_and_nutrients"
}

# Define the threshold for what is considered a "huge" and invalid yield value
HUGE_THRESHOLD = 1e11  # This is equivalent to 10e+10

# --- Function to Process a Directory ---
def process_directory(dir_name, dir_path):
    # Find all NetCDF files in the specified directory
    nc_files = glob.glob(os.path.join(dir_path, "wofost_results_*.nc"))

    if not nc_files:
        print(f"Error: No NetCDF files found in the directory: {dir_path}")
        return None

    # Initialize a list to store the results from each file
    results = []

    print(f"\nProcessing files in {dir_name}...")
    # Loop through each file found
    for file in nc_files:
        try:
            # --- 1. Extract Metadata from Filename ---
            filename = os.path.basename(file)
            parts = filename.split('_')
            # Assumes filename format: wofost_results_SCENARIO_sowing_DATE.nc
            scenario = parts[2]
            tsm1 = scenario[:2]
            tsm2 = scenario[2:]
            sow_date = parts[4].replace('.nc', '')

            # --- 2. Load Data and Prevent Overflow ---
            # Open the NetCDF file, disabling automatic time decoding to avoid warnings
            with xr.open_dataset(file, decode_times=False) as ds:
                # Extract the 'Yield_Average' variable and immediately cast it to float64.
                yield_avg = ds['Yield_Average'].values.astype(np.float64)

            # --- 3. Filter Out Invalid Data ---
            invalid_mask = (
                np.isnan(yield_avg) |
                np.isinf(yield_avg) |
                (yield_avg <= 0) |
                (yield_avg == -9999) |
                (yield_avg > HUGE_THRESHOLD)
            )

            # For detailed reporting, we can also count how many were huge
            huge_count = np.sum((yield_avg > HUGE_THRESHOLD) & ~np.isnan(yield_avg))

            # Apply the mask to get an array of only the valid yield values
            valid_yield = yield_avg[~invalid_mask]

            # --- 4. Calculate Mean and Store Results ---
            if valid_yield.size > 0:
                mean_yield = np.mean(valid_yield)
            else:
                mean_yield = np.nan

            # Append the results for this file to our list
            results.append({
                'tsm1': tsm1,
                'tsm2': tsm2,
                'sowing_date': sow_date,
                'mean_yield': mean_yield,
                'invalid_points': np.sum(invalid_mask),
                'huge_points': huge_count
            })

        except Exception as e:
            print(f"Could not process file {file}. Error: {e}")

    # Convert the list of results into a pandas DataFrame for easy analysis
    df = pd.DataFrame(results)

    # Sort the DataFrame for better readability
    df = df.sort_values(['tsm1', 'tsm2', 'sowing_date'])

    # Print the main results table
    print(f"\n--- Results Table for {dir_name} (Mean Yields with Invalid Points Filtered) ---")
    print(df.to_string())

    # Print a summary of invalid points for each run
    print(f"\n--- Summary of Invalid Grid Points for {dir_name} ---")
    for _, row in df.iterrows():
        print(f"Scenario {row['tsm1']}{row['tsm2']} with sow_date '{row['sowing_date']}': {row['invalid_points']} invalid points (of which {row['huge_points']} were > {HUGE_THRESHOLD:.0e})")

    return df

# --- Process All Directories ---
dfs = {}
for dir_name, dir_path in data_dirs.items():
    dfs[dir_name] = process_directory(dir_name, dir_path)

if any(df is None for df in dfs.values()):
    print("Error processing one or more directories. Exiting.")
    exit()

df_original = dfs["original"]
df_evtra = dfs["evtra"]
df_nutrients = dfs["nutrients"]
df_evtra_nutrients = dfs["evtra_nutrients"]

# --- Generate Plots for Each Directory ---
# Get a sorted list of unique sowing dates (assuming same across dirs)
unique_sows = sorted(df_original['sowing_date'].unique())

# Define the desired order and new labels for the axes
axis_order = ['e1', 'a1', 'l1']
axis_labels = ['early', 'average', 'late']

# For y-axis in heatmap: late at top, early at bottom
y_axis_order = ['l1', 'a1', 'e1']  # late (top), average, early (bottom)
y_axis_labels = ['late', 'average', 'early']

# Remap tsm1 and tsm2 to descriptive labels for plotting
label_map = {'e1': 'early', 'a1': 'average', 'l1': 'late'}

for dir_name, df in dfs.items():
    df['flowering'] = df['tsm1'].map(label_map)
    df['maturity'] = df['tsm2'].map(label_map)

    print(f"\nGenerating heatmaps for {dir_name}...")
    for sow in unique_sows:
        sub_df = df[df['sowing_date'] == sow]
        if not sub_df.empty:
            # Pivot the data and reorder the index and columns
            pivot = sub_df.pivot_table(index='tsm1', columns='tsm2', values='mean_yield')
            pivot = pivot.reindex(index=y_axis_order, columns=axis_order)

            plt.figure(figsize=(9, 8))

            # Create the heatmap
            ax = sns.heatmap(pivot, annot=True, fmt=".0f", cmap="RdYlGn", linewidths=.5,
                              cbar_kws={'label': 'Mean Yield (kg/ha)'})

            # Set titles and labels
            ax.set_title(f"Yield Heatmap for Sowing Date: '{sow}' ({dir_name})", pad=20)
            ax.set_xlabel('Maturity timing (parameter Tsum2)', labelpad=15)
            ax.set_ylabel('Flowering timing (parameter Tsum1)', labelpad=15)

            # Set the custom tick labels
            ax.set_xticklabels(axis_labels)
            ax.set_yticklabels(y_axis_labels, rotation=0)

            plt.tight_layout()
            plt.savefig(f'yield_heatmap_{dir_name}_sow_{sow}.png')
            print(f"Heatmap saved as 'yield_heatmap_{dir_name}_sow_{sow}.png'")
            plt.close()

    print(f"\nGenerating contour maps for {dir_name}...")
    for sow in unique_sows:
        sub_df = df[df['sowing_date'] == sow]
        pivot = sub_df.pivot_table(index='flowering', columns='maturity', values='mean_yield')
        # Order for contour consistent with heatmap: late top, early bottom - use same order but invert y-axis
        pivot = pivot.reindex(index=['late', 'average', 'early'], columns=['early', 'average', 'late'])

        plt.figure(figsize=(8, 6))
        sns.set_theme(style="white")
        levels = np.linspace(pivot.min().min(), pivot.max().max(), 10)
        contour = plt.contourf(pivot.columns, pivot.index, pivot.values, levels=levels, cmap='RdYlGn')
        plt.colorbar(contour, label='Mean Yield (kg/ha)')
        plt.title(f"Yield Contour for Sowing: {sow} ({dir_name})", fontsize=14, pad=20)
        plt.xlabel('Maturity Timing', fontsize=12)
        plt.ylabel('Flowering Timing', fontsize=12, labelpad=20)
        plt.grid(False)
        # Invert y-axis to have 'late' at top, 'early' at bottom
        ax = plt.gca()
        ax.invert_yaxis()
        ax.tick_params(axis='y', pad=30)  # Move tick labels more to the left
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                value = pivot.values[i, j]
                plt.text(pivot.columns[j], pivot.index[i], f"{value:.0f}", ha='center', va='center', 
                         color='black', fontsize=14, fontweight='bold',
                         path_effects=[pe.withStroke(linewidth=1.5, foreground='white')])
        plt.tight_layout()
        plt.subplots_adjust(left=0.22)  # Adjusted space on left
        plt.savefig(f'yield_contour_{dir_name}_sow_{sow}.png', dpi=300)
        print(f"Contour plot saved as 'yield_contour_{dir_name}_sow_{sow}.png'")
        plt.close()

# --- Comparison Between Directories ---
print("\n--- Comparison Between Evtra and Original ---")

# Merge the two dataframes on common keys
df_comparison_evtra = df_original[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']].merge(
    df_evtra[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']],
    on=['tsm1', 'tsm2', 'sowing_date'],
    suffixes=('_original', '_evtra')
)
df_comparison_evtra['difference'] = df_comparison_evtra['mean_yield_evtra'] - df_comparison_evtra['mean_yield_original']

# Print comparison table
print("\n--- Yield Differences (Evtra - Original) ---")
print(df_comparison_evtra.to_string())

# Create text-based pivot tables for differences
for sow in unique_sows:
    sub_df = df_comparison_evtra[df_comparison_evtra['sowing_date'] == sow]
    if not sub_df.empty:
        pivot_diff = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference').round(1)
        print(f"\n--- Yield Difference (kg/ha) for sowing_date '{sow}' (Evtra - Original) ---")
        print("(rows: tsm1, columns: tsm2)")
        print(pivot_diff)

# Generate difference heatmaps
print("\nGenerating difference heatmaps for Evtra vs Original...")
for sow in unique_sows:
    sub_df = df_comparison_evtra[df_comparison_evtra['sowing_date'] == sow]
    if not sub_df.empty:
        pivot = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference')
        pivot = pivot.reindex(index=y_axis_order, columns=axis_order)

        plt.figure(figsize=(9, 8))

        # Create the heatmap for differences (use diverging colormap)
        ax = sns.heatmap(pivot, annot=True, fmt=".0f", cmap="RdBu", linewidths=.5, center=0,
                          cbar_kws={'label': 'Yield Difference (kg/ha)'})

        # Set titles and labels
        ax.set_title(f"Yield Difference Heatmap for Sowing Date: '{sow}' (Evtra - Original)", pad=20)
        ax.set_xlabel('Maturity timing (parameter Tsum2)', labelpad=15)
        ax.set_ylabel('Flowering timing (parameter Tsum1)', labelpad=15)

        # Set the custom tick labels
        ax.set_xticklabels(axis_labels)
        ax.set_yticklabels(y_axis_labels, rotation=0)

        plt.tight_layout()
        plt.savefig(f'yield_diff_heatmap_evtra_sow_{sow}.png')
        print(f"Difference heatmap saved as 'yield_diff_heatmap_evtra_sow_{sow}.png'")
        plt.close()

# Comparison for Nutrients vs Original
print("\n--- Comparison Between Nutrients and Original ---")

df_comparison_nutrients = df_original[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']].merge(
    df_nutrients[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']],
    on=['tsm1', 'tsm2', 'sowing_date'],
    suffixes=('_original', '_nutrients')
)
df_comparison_nutrients['difference'] = df_comparison_nutrients['mean_yield_nutrients'] - df_comparison_nutrients['mean_yield_original']

# Print comparison table
print("\n--- Yield Differences (Nutrients - Original) ---")
print(df_comparison_nutrients.to_string())

# Create text-based pivot tables for differences
for sow in unique_sows:
    sub_df = df_comparison_nutrients[df_comparison_nutrients['sowing_date'] == sow]
    if not sub_df.empty:
        pivot_diff = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference').round(1)
        print(f"\n--- Yield Difference (kg/ha) for sowing_date '{sow}' (Nutrients - Original) ---")
        print("(rows: tsm1, columns: tsm2)")
        print(pivot_diff)

# Generate difference heatmaps
print("\nGenerating difference heatmaps for Nutrients vs Original...")
for sow in unique_sows:
    sub_df = df_comparison_nutrients[df_comparison_nutrients['sowing_date'] == sow]
    if not sub_df.empty:
        pivot = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference')
        pivot = pivot.reindex(index=y_axis_order, columns=axis_order)

        plt.figure(figsize=(9, 8))

        # Create the heatmap for differences (use diverging colormap)
        ax = sns.heatmap(pivot, annot=True, fmt=".0f", cmap="RdBu", linewidths=.5, center=0,
                          cbar_kws={'label': 'Yield Difference (kg/ha)'})

        # Set titles and labels
        ax.set_title(f"Yield Difference Heatmap for Sowing Date: '{sow}' (Nutrients - Original)", pad=20)
        ax.set_xlabel('Maturity timing (parameter Tsum2)', labelpad=15)
        ax.set_ylabel('Flowering timing (parameter Tsum1)', labelpad=15)

        # Set the custom tick labels
        ax.set_xticklabels(axis_labels)
        ax.set_yticklabels(y_axis_labels, rotation=0)

        plt.tight_layout()
        plt.savefig(f'yield_diff_heatmap_nutrients_sow_{sow}.png')
        print(f"Difference heatmap saved as 'yield_diff_heatmap_nutrients_sow_{sow}.png'")
        plt.close()

# Comparison for Evtra_Nutrients vs Original
print("\n--- Comparison Between Evtra_Nutrients and Original ---")

df_comparison_evtra_nutrients = df_original[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']].merge(
    df_evtra_nutrients[['tsm1', 'tsm2', 'sowing_date', 'mean_yield']],
    on=['tsm1', 'tsm2', 'sowing_date'],
    suffixes=('_original', '_evtra_nutrients')
)
df_comparison_evtra_nutrients['difference'] = df_comparison_evtra_nutrients['mean_yield_evtra_nutrients'] - df_comparison_evtra_nutrients['mean_yield_original']

# Print comparison table
print("\n--- Yield Differences (Evtra_Nutrients - Original) ---")
print(df_comparison_evtra_nutrients.to_string())

# Create text-based pivot tables for differences
for sow in unique_sows:
    sub_df = df_comparison_evtra_nutrients[df_comparison_evtra_nutrients['sowing_date'] == sow]
    if not sub_df.empty:
        pivot_diff = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference').round(1)
        print(f"\n--- Yield Difference (kg/ha) for sowing_date '{sow}' (Evtra_Nutrients - Original) ---")
        print("(rows: tsm1, columns: tsm2)")
        print(pivot_diff)

# Generate difference heatmaps
print("\nGenerating difference heatmaps for Evtra_Nutrients vs Original...")
for sow in unique_sows:
    sub_df = df_comparison_evtra_nutrients[df_comparison_evtra_nutrients['sowing_date'] == sow]
    if not sub_df.empty:
        pivot = sub_df.pivot_table(index='tsm1', columns='tsm2', values='difference')
        pivot = pivot.reindex(index=y_axis_order, columns=axis_order)

        plt.figure(figsize=(9, 8))

        # Create the heatmap for differences (use diverging colormap)
        ax = sns.heatmap(pivot, annot=True, fmt=".0f", cmap="RdBu", linewidths=.5, center=0,
                          cbar_kws={'label': 'Yield Difference (kg/ha)'})

        # Set titles and labels
        ax.set_title(f"Yield Difference Heatmap for Sowing Date: '{sow}' (Evtra_Nutrients - Original)", pad=20)
        ax.set_xlabel('Maturity timing (parameter Tsum2)', labelpad=15)
        ax.set_ylabel('Flowering timing (parameter Tsum1)', labelpad=15)

        # Set the custom tick labels
        ax.set_xticklabels(axis_labels)
        ax.set_yticklabels(y_axis_labels, rotation=0)

        plt.tight_layout()
        plt.savefig(f'yield_diff_heatmap_evtra_nutrients_sow_{sow}.png')
        print(f"Difference heatmap saved as 'yield_diff_heatmap_evtra_nutrients_sow_{sow}.png'")
        plt.close()

print("\n--- Analysis Complete ---")
