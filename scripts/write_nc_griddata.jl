using NCDatasets, Dates, ProgressMeter, LinearAlgebra

# --- Configuration ---
const CROP_MASK_PATH = "./springwheat-char-05d_C3S-glob-agric_2005_v5.nc"
const TAIR_PATH = "./Tair_1989-2019.nc"
const OUTPUT_PATH = "./all_griddata.nc"
const START_YEAR = 1989
const END_YEAR = 2019

# Define all 9 combinations of sowing and maturation variables
const COMBINATIONS = [
    (sow="sow_e1", mat="mat_e1", suffix="e1e1"),
    (sow="sow_e1", mat="mat_a1", suffix="e1a1"),
    (sow="sow_e1", mat="mat_l1", suffix="e1l1"),
    (sow="sow_a1", mat="mat_e1", suffix="a1e1"),
    (sow="sow_a1", mat="mat_a1", suffix="a1a1"),
    (sow="sow_a1", mat="mat_l1", suffix="a1l1"),
    (sow="sow_l1", mat="mat_e1", suffix="l1e1"),
    (sow="sow_l1", mat="mat_a1", suffix="l1a1"),
    (sow="sow_l1", mat="mat_l1", suffix="l1l1"),
]

# Dekad mapping
const DEKADS = (
    (m=1,d=10), (m=1,d=20), (m=1,d=31), (m=2,d=10), (m=2,d=20), (m=2,d=28),
    (m=3,d=10), (m=3,d=20), (m=3,d=31), (m=4,d=10), (m=4,d=20), (m=4,d=30),
    (m=5,d=10), (m=5,d=20), (m=5,d=31), (m=6,d=10), (m=6,d=20), (m=6,d=30),
    (m=7,d=10), (m=7,d=20), (m=7,d=31), (m=8,d=10), (m=8,d=20), (m=8,d=30),
    (m=9,d=10), (m=9,d=20), (m=9,d=30), (m=10,d=10), (m=10,d=20), (m=10,d=31),
    (m=11,d=10), (m=11,d=20), (m=11,d=30), (m=12,d=10), (m=12,d=20), (m=12,d=31)
)

function main()
    println("Loading input files...")

    # --- Load Data ---
    ds_mask = NCDataset(CROP_MASK_PATH)
    ds_tair = NCDataset(TAIR_PATH)

    # --- FIX: Removed transpose() ---
    # In Julia, NCDatasets loads dimensions as (lon, lat, time)
    sow_vars = Dict(v => ds_mask[v][:,:,1] for v in ["sow_e1", "sow_a1", "sow_l1"])
    mat_vars = Dict(v => ds_mask[v][:,:,1] for v in ["mat_e1", "mat_a1", "mat_l1"])
    tsum_ea = ds_mask["tsumEA"][:,:,1]
    tsum_am = ds_mask["tsumAM"][:,:,1]
    
    # Immediately copy attributes into a new Dictionary while file is open
    sow_attribs = Dict(v => Dict(k => ds_mask[v].attrib[k] for k in keys(ds_mask[v].attrib)) for v in ["sow_e1", "sow_a1", "sow_l1"])
    
    lat = ds_mask["lat"][:]
    lon = ds_mask["lon"][:]

    temp_celsius = ds_tair["Tair"][:,:,:] .- 273.15
    dates = ds_tair["time"][:]

    close(ds_mask)
    close(ds_tair)
    println("Loading complete...")

    # --- Create NetCDF file ---
    NCDataset(OUTPUT_PATH, "c") do ds
        println("Defining NetCDF dimensions...")
        defDim(ds, "lon", length(lon))
        defDim(ds, "lat", length(lat))
        defVar(ds, "lon", lon, ("lon",))
        defVar(ds, "lat", lat, ("lat",))

        # Copy sow variables and their metadata
        println("Copying sowing variables from source...")
        for v_name in ["sow_e1", "sow_a1", "sow_l1"]
            attrib_dict = sow_attribs[v_name]
            fill_val = pop!(attrib_dict, "_FillValue", -9999.0f0)
            defVar(ds, v_name, sow_vars[v_name], ("lon", "lat"), attrib=attrib_dict, fillvalue=fill_val)
        end

        # --- Main Calculation and Writing Loop ---
        p_outer = Progress(length(COMBINATIONS), 1, "Overall Progress: ")

        for combo in COMBINATIONS
            ProgressMeter.update!(p_outer, desc="Processing combination: $(combo.suffix)...")
            sow_data = sow_vars[combo.sow]
            mat_data = mat_vars[combo.mat]

            tsum_grid = fill(NaN, length(lon), length(lat))
            tsum1_grid = fill(NaN, length(lon), length(lat))
            tsum2_grid = fill(NaN, length(lon), length(lat))

            p_inner = Progress(length(lon) * length(lat), 0.1, "  ↳ Analyzing cells for $(combo.suffix): ")

            for ln in 1:length(lon), lt in 1:length(lat)
                # --- FIX: Changed indexing from [lt, ln] to [ln, lt] ---
                sow_dekad = sow_data[ln, lt]
                mat_dekad = mat_data[ln, lt]

                if !ismissing(sow_dekad) && sow_dekad > 0 && !ismissing(mat_dekad)
                    sow_dekad_int = round(Int, sow_dekad)
                    mat_duration_days = 10 * mat_dekad
                    total_tsum = 0.0
                    year_count = 0

                    for year_val in START_YEAR:END_YEAR
                        dekad_data = DEKADS[sow_dekad_int]
                        emergence = Date(year_val, dekad_data.m, dekad_data.d)
                        maturity = emergence + Day(mat_duration_days)
                        season_indices = findall(d -> emergence <= Date(d) <= maturity, dates)

                        if !isempty(season_indices) && Date(maturity) <= dates[end]
                            yearly_tsum = sum(skipmissing(view(temp_celsius, ln, lt, season_indices)))
                            total_tsum += yearly_tsum
                            year_count += 1
                        end
                    end

                    if year_count >= 30 && total_tsum > 100.0
                        avg_tsum = round(Int, total_tsum / year_count)
                        tsum_grid[ln, lt] = avg_tsum
                        # --- FIX: Changed indexing from [lt, ln] to [ln, lt] ---
                        cell_tsum_ea = tsum_ea[ln, lt]
                        cell_tsum_am = tsum_am[ln, lt]

                        if !ismissing(cell_tsum_ea) && !ismissing(cell_tsum_am)
                            total_pheno_tsum = cell_tsum_ea + cell_tsum_am
                            if total_pheno_tsum > 0
                                tsum1_grid[ln, lt] = round(Int, (cell_tsum_ea / total_pheno_tsum) * avg_tsum)
                                tsum2_grid[ln, lt] = round(Int, (cell_tsum_am / total_pheno_tsum) * avg_tsum)
                            end
                        end
                    end
                end
                ProgressMeter.next!(p_inner)
            end

            defVar(ds, "avg_tsum_$(combo.suffix)", tsum_grid, ("lon", "lat"); fillvalue=NaN)
            defVar(ds, "avg_tsum1_$(combo.suffix)", tsum1_grid, ("lon", "lat"); fillvalue=NaN)
            defVar(ds, "avg_tsum2_$(combo.suffix)", tsum2_grid, ("lon", "lat"); fillvalue=NaN)

            ProgressMeter.next!(p_outer)
        end
    end

    println("\nProcessing complete. Output written to $OUTPUT_PATH")
end

main()
