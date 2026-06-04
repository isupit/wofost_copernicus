#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "wofost.h"
#include "extern.h"
#include "output_netcdf.h"
#include "input_griddata.h"
#include "input_fertilizer.h"

// --- Global variables for optional flags ---
// Initialize them to their default state (0 = off/false)
int use_potential_nutrients = 0;
int use_potential_evtra = 0;
int use_gridded_tsum = 0;
float tsum1_offset = 0.0f;
float tsum2_offset = 0.0f;

int main(int argc, char **argv)
{
    FILE **files; 
    FILE *fptr;   

    SimUnit *initial = NULL; 
    Weather *head;          
    Green *wipe;             

    int CycleLength = 300; 
    int NumberOfFiles;     
    int Emergence;         
    int i;               

    char list[MAX_STRING];      
    char meteolist[MAX_STRING]; 
    char name[MAX_STRING];      
    char name_old[MAX_STRING]; 

    /* --- Variable names for gridded crop data input --- */
    char grid_data_file[MAX_STRING];
    char tsum1_var[MAX_STRING]; 
    char tsum2_var[MAX_STRING]; 
    char sow_var[MAX_STRING];  
    
    char n_fert_file[MAX_STRING];

    char output_file[MAX_STRING]; /* Dynamic output filename */

    Step = 1.; 

    // Initialize all strings to empty
    memset(list, '\0', MAX_STRING);
    memset(meteolist, '\0', MAX_STRING);
    memset(grid_data_file, '\0', MAX_STRING); 
    memset(tsum1_var, '\0', MAX_STRING); 
    memset(tsum2_var, '\0', MAX_STRING); 
    memset(sow_var, '\0', MAX_STRING);
    memset(n_fert_file, '\0', MAX_STRING);
    
    // IMPROVED: Flexible argument parsing that handles all modes
    if (argc < 3) {
        goto usage_error;
    }
    
    // Always get the first two mandatory arguments
    if (strlen(argv[1]) >= MAX_STRING || strlen(argv[2]) >= MAX_STRING) {
        fprintf(stderr, "Error: Argument too long\n");
        exit(1);
    }
    strncpy(list, argv[1], strlen(argv[1]));
    strncpy(meteolist, argv[2], strlen(argv[2]));
    
    // Parse remaining arguments based on content
    int arg_idx = 3;
    int tsum_provided = 0;
    int sow_provided = 0;
    int grid_file_provided = 0;
    
    use_gridded_tsum = 0;
    
    // Parse positional arguments until we hit flags
    while (arg_idx < argc && argv[arg_idx][0] != '-') {
        if (!grid_file_provided) {
            // First extra argument is the grid data file (if it looks like a filename)
            if (strstr(argv[arg_idx], ".nc") != NULL || strstr(argv[arg_idx], ".nc4") != NULL) {
                strncpy(grid_data_file, argv[arg_idx], strlen(argv[arg_idx]));
                grid_file_provided = 1;
            } else {
                // Not a grid file, treat as sowing variable (mixed mode)
                strncpy(sow_var, argv[arg_idx], strlen(argv[arg_idx]));
                sow_provided = 1;
            }
        } else if (strstr(argv[arg_idx], "tsum") != NULL) {
            // Arguments containing "tsum" are TSUM variables
            if (tsum_provided == 0) {
                strncpy(tsum1_var, argv[arg_idx], strlen(argv[arg_idx]));
                tsum_provided++;
            } else if (tsum_provided == 1) {
                strncpy(tsum2_var, argv[arg_idx], strlen(argv[arg_idx]));
                tsum_provided++;
            }
        } else if (strstr(argv[arg_idx], "sow") != NULL) {
            // Arguments containing "sow" are sowing variables
            strncpy(sow_var, argv[arg_idx], strlen(argv[arg_idx]));
            sow_provided = 1;
        } else {
            // Unknown positional argument
            fprintf(stderr, "Error: Unrecognized positional argument '%s'\n", argv[arg_idx]);
            goto usage_error;
        }
        arg_idx++;
    }
    
    // Determine mode based on what we found
    if (tsum_provided == 2 && strlen(sow_var) > 0 && grid_file_provided) {
        // Full gridded TSUM mode
        use_gridded_tsum = 1;
        printf("Full gridded mode: TSUM1=%s, TSUM2=%s, SOW=%s from %s\n", tsum1_var, tsum2_var, sow_var, grid_data_file);
    } else if (strlen(sow_var) > 0 && grid_file_provided) {
        // Mixed mode: default TSUM + gridded sowing
        use_gridded_tsum = 0;
        printf("Mixed mode: default TSUM + gridded sowing (%s) from %s\n", sow_var, grid_data_file);
    } else if (strlen(sow_var) > 0) {
        // Sowing variable name only (no grid file) - use default sowing dates
        use_gridded_tsum = 0;
        printf("Default TSUM mode: sowing variable '%s' specified but no grid file - using default sowing dates\n", sow_var);
    } else if (grid_file_provided) {
        // Only grid file provided - error
        fprintf(stderr, "Error: Grid file '%s' provided but no TSUM or sowing variables specified\n", grid_data_file);
        goto usage_error;
    } else {
        // No grid components - pure default mode
        use_gridded_tsum = 0;
        printf("Default mode: using crop parameter defaults for TSUM and sowing\n");
    }
    
    // Parse remaining arguments as flags (starting from current arg_idx)
    for (int i = arg_idx; i < argc; i++) {
        if (strcmp(argv[i], "--use-potential-nutrients") == 0) {
            use_potential_nutrients = 1;
        } else if (strcmp(argv[i], "--use-potential-evtra") == 0) {
            use_potential_evtra = 1;
        } else if (strcmp(argv[i], "--n-fertilizer-nc") == 0) {
            if (i + 1 < argc) {
                strncpy(n_fert_file, argv[i + 1], MAX_STRING - 1);
                i++; // Skip the filename
            } else {
                fprintf(stderr, "Error: --n-fertilizer-nc flag requires a filename.\n");
                exit(1);
            }
        } else if (strcmp(argv[i], "--tsum1-offset") == 0) {
            if (i + 1 < argc) {
                tsum1_offset = atof(argv[i + 1]);
                i++; // Skip the value
            } else {
                fprintf(stderr, "Error: --tsum1-offset flag requires a float value.\n");
                exit(1);
            }
        } else if (strcmp(argv[i], "--tsum2-offset") == 0) {
            if (i + 1 < argc) {
                tsum2_offset = atof(argv[i + 1]);
                i++; // Skip the value
            } else {
                fprintf(stderr, "Error: --tsum2-offset flag requires a float value.\n");
                exit(1);
            }
        } else {
            fprintf(stderr, "Error: Unknown flag '%s'\n", argv[i]);
            goto usage_error;
        }
    }

    // Configuration summary
    printf("\n--- Configuration Summary ---\n");
    printf("Simulation list: %s\n", list);
    printf("Meteo list: %s\n", meteolist);
    
    if (use_gridded_tsum) {
        printf("MODE: Gridded TSUM values\n");
        printf("  Grid data file: %s\n", grid_data_file);
        printf("  TSUM1 variable: %s\n", tsum1_var);
        printf("  TSUM2 variable: %s\n", tsum2_var);
        printf("  Sowing variable: %s\n", sow_var);
    } else {
        printf("MODE: Default TSUM values from crop parameter file\n");
    }
    
    printf("Optional flag --use-potential-nutrients set: %s\n", use_potential_nutrients ? "Yes" : "No");
    printf("Optional flag --use-potential-evtra set: %s\n", use_potential_evtra ? "Yes" : "No");
    printf("Optional N Fertilizer NetCDF provided: %s\n", strlen(n_fert_file) > 0 ? n_fert_file : "No");
    printf("TSUM1 offset: %.1f\n", tsum1_offset);
    printf("TSUM2 offset: %.1f\n", tsum2_offset);
    printf("---------------------------\n\n");

    /* --- Construct dynamic output filename --- */
    char offset_suffix[128] = "";
    char base_suffix[128] = "";
    
    // Build the offset suffix with clear labels - ALWAYS INCLUDE VALUES
    char t1_part[64] = "";
    char t2_part[64] = "";
    
    // Build t1 part (ALWAYS include, even if zero)
    {
        char tsum1_str[16];
        snprintf(tsum1_str, sizeof(tsum1_str), "%d", (int)tsum1_offset);
        if (tsum1_offset > 0) {
            snprintf(t1_part, sizeof(t1_part), "_t1+%s", tsum1_str);
        } else if (tsum1_offset < 0) {
            snprintf(t1_part, sizeof(t1_part), "_t1%s", tsum1_str);  // Negative sign already in number
        } else {
            snprintf(t1_part, sizeof(t1_part), "_t1+0");
        }
    }
    
    // Build t2 part (ALWAYS include, even if zero)  
    {
        char tsum2_str[16];
        snprintf(tsum2_str, sizeof(tsum2_str), "%d", (int)tsum2_offset);
        if (tsum2_offset > 0) {
            snprintf(t2_part, sizeof(t2_part), "_t2+%s", tsum2_str);
        } else if (tsum2_offset < 0) {
            snprintf(t2_part, sizeof(t2_part), "_t2%s", tsum2_str);  // Negative sign already in number
        } else {
            snprintf(t2_part, sizeof(t2_part), "_t2+0");
        }
    }
    
    // ALWAYS combine t1 and t2 parts into offset_suffix
    snprintf(offset_suffix, sizeof(offset_suffix), "%s%s", t1_part, t2_part);
    
    // Now build the base filename
    if (use_gridded_tsum) {
        /* Full gridded TSUM mode */
        if (strlen(tsum1_var) == 0) {
            fprintf(stderr, "Error: tsum1_var is empty when using gridded TSUM\n");
            exit(1);
        }
        if (strlen(sow_var) == 0) {
            fprintf(stderr, "Error: sow_var is empty when using gridded TSUM\n");
            exit(1);
        }
        
        char *tsum_suffix = strrchr(tsum1_var, '_');
        if (tsum_suffix && strlen(tsum_suffix) > 1) {
            tsum_suffix++; // Skip the '_'
            snprintf(base_suffix, sizeof(base_suffix), "%s_%s", tsum_suffix, sow_var);
        } else {
            snprintf(base_suffix, sizeof(base_suffix), "default_%s", sow_var);
        }
        
    } else if (strlen(sow_var) > 0 && strlen(grid_data_file) > 0) {
        /* Mixed mode: default TSUM + gridded sowing */
        char *sow_suffix = strrchr(sow_var, '_');
        if (sow_suffix && strlen(sow_suffix) > 1) {
            sow_suffix++; // Skip the '_'
            snprintf(base_suffix, sizeof(base_suffix), "default_tsum_sow_%s", sow_suffix);
        } else {
            snprintf(base_suffix, sizeof(base_suffix), "default_tsum_sow_%s", sow_var);
        }
        
    } else {
        /* Pure default mode */
        snprintf(base_suffix, sizeof(base_suffix), "default");
    }
    
    /* --- Create per-run output directory: output/YYYYMMDD_HHMMSS_<suffix>/ --- */
    time_t now = time(NULL);
    struct tm *lt = localtime(&now);
    char run_timestamp[32];
    strftime(run_timestamp, sizeof(run_timestamp), "%Y%m%d_%H%M%S", lt);

    char run_dir[MAX_STRING];
    snprintf(run_dir, MAX_STRING, "output/%s_%s%s", run_timestamp, base_suffix, offset_suffix);

    /* Create output/ then the run subdir (ignore errors if they already exist) */
    mkdir("output", 0755);
    if (mkdir(run_dir, 0755) != 0) {
        fprintf(stderr, "Warning: could not create run directory '%s'\n", run_dir);
    }

    /* Full path for the NetCDF output file */
    snprintf(output_file, MAX_STRING, "%s/wofost_results_%s%s.nc", run_dir, base_suffix, offset_suffix);

    printf("Constructed output file: %s\n", output_file);

    /* --- Write run.log with all inputs for reproducibility --- */
    {
        char logpath[MAX_STRING];
        snprintf(logpath, MAX_STRING, "%s/run.log", run_dir);
        FILE *logfp = fopen(logpath, "w");
        if (logfp == NULL) {
            fprintf(stderr, "Warning: could not create run log '%s'\n", logpath);
        } else {
            /* Timestamp */
            char timebuf[64];
            strftime(timebuf, sizeof(timebuf), "%Y-%m-%d %H:%M:%S", lt);
            fprintf(logfp, "=== WOFOST Run Log ===\n");
            fprintf(logfp, "Timestamp       : %s\n", timebuf);
#ifdef GIT_HASH
            fprintf(logfp, "Git commit      : %s\n", GIT_HASH);
#else
            fprintf(logfp, "Git commit      : (unknown - not compiled with git hash)\n");
#endif
            fprintf(logfp, "\n");

            /* Reconstruct command line */
            fprintf(logfp, "=== Command Line ===\n");
            for (int ci = 0; ci < argc; ci++) {
                fprintf(logfp, "%s%s", argv[ci], ci < argc - 1 ? " " : "\n");
            }
            fprintf(logfp, "\n");

            /* Configuration summary */
            fprintf(logfp, "=== Configuration ===\n");
            fprintf(logfp, "Simulation list         : %s\n", list);
            fprintf(logfp, "Meteo list              : %s\n", meteolist);
            if (use_gridded_tsum) {
                fprintf(logfp, "Mode                    : Full gridded TSUM\n");
                fprintf(logfp, "Grid data file          : %s\n", grid_data_file);
                fprintf(logfp, "TSUM1 variable          : %s\n", tsum1_var);
                fprintf(logfp, "TSUM2 variable          : %s\n", tsum2_var);
                fprintf(logfp, "Sowing variable         : %s\n", sow_var);
            } else if (strlen(sow_var) > 0 && strlen(grid_data_file) > 0) {
                fprintf(logfp, "Mode                    : Default TSUM + gridded sowing\n");
                fprintf(logfp, "Grid data file          : %s\n", grid_data_file);
                fprintf(logfp, "Sowing variable         : %s\n", sow_var);
            } else {
                fprintf(logfp, "Mode                    : Default TSUM + default sowing\n");
            }
            fprintf(logfp, "TSUM1 offset            : %.1f\n", tsum1_offset);
            fprintf(logfp, "TSUM2 offset            : %.1f\n", tsum2_offset);
            fprintf(logfp, "Use potential nutrients  : %s\n", use_potential_nutrients ? "Yes" : "No");
            fprintf(logfp, "Use potential evtra     : %s\n", use_potential_evtra ? "Yes" : "No");
            fprintf(logfp, "N fertilizer file       : %s\n", strlen(n_fert_file) > 0 ? n_fert_file : "(none)");
            fprintf(logfp, "Output directory        : %s\n", run_dir);
            fprintf(logfp, "Output NetCDF           : %s\n", output_file);
            fprintf(logfp, "\n");

            /* Dump list.txt */
            fprintf(logfp, "=== Simulation List (%s) ===\n", list);
            FILE *lf = fopen(list, "r");
            if (lf) {
                char lbuf[MAX_STRING];
                while (fgets(lbuf, MAX_STRING, lf)) fprintf(logfp, "%s", lbuf);
                fclose(lf);
            } else {
                fprintf(logfp, "(could not open)\n");
            }
            fprintf(logfp, "\n");

            /* Dump meteolist.txt */
            fprintf(logfp, "=== Meteo List (%s) ===\n", meteolist);
            FILE *mf = fopen(meteolist, "r");
            if (mf) {
                char mbuf[MAX_STRING];
                while (fgets(mbuf, MAX_STRING, mf)) fprintf(logfp, "%s", mbuf);
                fclose(mf);
            } else {
                fprintf(logfp, "(could not open)\n");
            }
            fprintf(logfp, "\n");

            fclose(logfp);
            printf("Run log written to: %s\n", logpath);
        }
    }

    /* Fill the crop, soil, site and management place holders*/ 
    NumberOfFiles = GetSimInput(list);

    /* Set the initial Grid address */ 
    initial = Grid;

    /* Get the meteo filenames and put them in the placeholder */ 
    GetMeteoInput(meteolist);

    /* Allocate memory for the file pointers */      
    files = malloc(sizeof(**files) * NumberOfFiles); 

    /* Go back to the beginning of the list */ 
    Grid = initial;

    /* Open the output files */                            
    memset(name_old, '\0', MAX_STRING);                    
    while (Grid)                                           
    {                                                      
        memset(name, '\0', MAX_STRING);                    
        strncpy(name, Grid->output, strlen(Grid->output)); 

        if (strcmp(name_old, name) != 0) 
        {
            files[Grid->file] = fptr = fopen(name, "w");
            if (files[Grid->file] == NULL) 
            {
                fprintf(stderr, "Cannot initialize output file %s.\n", name); 
                exit(0);                                              
            }
            header(files[Grid->file]); 
        }
        else
        {
            if (fptr != NULL)
                files[Grid->file] = fptr;
            else
            {
                fprintf(stderr, "Cannot initialize file pointer\n"); 
                exit(0);                                            
            }
        }

        // allocate memory for the statistical analysis 
        for (i = 0; i <= Meteo->Seasons; i++)
        {
            Grid->twso[i] = 0.0;
            Grid->length[i] = 0.0;
            Grid->applied_n[i]= 0.0f;
            Grid->cold_days[i] = 0;   /* NEW */
        }

        memset(name_old, '\0', MAX_STRING);                    
        strncpy(name_old, Grid->output, strlen(Grid->output)); 

        Grid = Grid->next; 
    }

    // Go back to the beginning of the list
    Grid = initial;

    /* --- Setup NetCDF file --- */
    NcFile nc_output;
    printf("Setting up NetCDF output file '%s'...\n", output_file);
    // Use sow_var even in mixed mode, or empty string in pure default mode
    char sow_for_nc[MAX_STRING];
    if (strlen(sow_var) > 0) {
        strncpy(sow_for_nc, sow_var, MAX_STRING - 1);
    } else {
        strcpy(sow_for_nc, "default");
    }

    while (Meteo)
    {
        /* Get the meteodata */ 
        if (GetMeteoData(Meteo) != 1)
        {
            fprintf(stderr, "Cannot get meteo data.\n");
            exit(0);
        }

        SetupNetCDF(output_file, &nc_output, Meteo->nlat, Meteo->nlon, Meteo->Seasons, tsum1_var, tsum2_var, sow_for_nc);

        /* Load crop grid data conditionally */
        if (use_gridded_tsum) {
            printf("Using gridded TSUM1/TSUM2 data from %s\n", grid_data_file);
            GetGridData(Meteo, grid_data_file, tsum1_var, tsum2_var, sow_var);
        } else {
            printf("Using default TSUM1/TSUM2 values from crop parameter file\n");
            SetDefaultTSUM(Meteo);
            
            // Load sowing dates ONLY if we have both a grid file AND sowing variable
            if (strlen(sow_var) > 0 && strlen(grid_data_file) > 0) {
                printf("Loading gridded sowing dates from %s (%s variable)\n", grid_data_file, sow_var);
                LoadSowingDateOnly(Meteo, grid_data_file, sow_var);
            } else if (strlen(sow_var) > 0) {
                // Just sowing variable name provided, no grid file - use default sowing
                printf("Sowing variable '%s' provided but no grid file - using default sowing date (April 1st)\n", sow_var);
                for (size_t i = 0; i < Meteo->nlat; i++) {
                    for (size_t j = 0; j < Meteo->nlon; j++) {
                        Meteo->sowing_date_grid[i][j] = 10.0f;  // April 1st (dekad 10)
                    }
                }
            } else {
                // No sowing info - use default
                printf("No sowing information provided - using default sowing date (April 1st)\n");
                for (size_t i = 0; i < Meteo->nlat; i++) {
                    for (size_t j = 0; j < Meteo->nlon; j++) {
                        Meteo->sowing_date_grid[i][j] = 10.0f;  // April 1st (dekad 10)
                    }
                }
            }
        }
    

        /* --- Apply TSUM offsets to all grid cells --- */
        if (tsum1_offset != 0.0f || tsum2_offset != 0.0f) {
            printf("Applying TSUM offsets: TSUM1 +%.1f, TSUM2 +%.1f\n", tsum1_offset, tsum2_offset);
            ApplyTSUMOffsets(Meteo);
        }
    
        /* --- Load fertilizer data if provided --- */
        if (strlen(n_fert_file) > 0) {
            if (GetFertilizerData(Meteo, n_fert_file, "Total_inorg_N_application_rate") != 1) {
                fprintf(stderr, "Could not load N fertilizer data.\n");
                exit(1);
            }
        }

        printf("running %d - %d\n", Meteo->StartYear, Meteo->EndYear);

        for (Lon = 0; Lon < Meteo->nlon; Lon++)
        {
            for (Lat = 0; Lat < Meteo->nlat; Lat++)
            {
                if (Mask[Lon][Lat] != 1) 
                {
                    continue;
                }
                
                /* Update Tsum values for this grid cell */
                Grid = initial;
                while (Grid) {
                    /* Update Tsum based on grid position (already includes offsets) */
                    Grid->crp->prm.TempSum1 = Meteo->tsum1_grid[Lat][Lon];
                    Grid->crp->prm.TempSum2 = Meteo->tsum2_grid[Lat][Lon];
                    
                    /* Convert day of year to "MM-DD" string */
                    int day_of_year = (int)Meteo->sowing_date_grid[Lat][Lon];
                    if (day_of_year < 1 || day_of_year > 365) {
                        strncpy(Grid->start, "01-01", 5);
                        Grid->start[5] = '\0';
                    } else {
                        int days_in_months[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
                        int month = 1;
                        int day = day_of_year;
                        for (int m = 0; m < 12; m++) {
                            if (day <= days_in_months[m]) {
                                month = m + 1;
                                break;
                            }
                            day -= days_in_months[m];
                        }
                        char date_str[6];
                        sprintf(date_str, "%02d-%02d", month, day);
                        strncpy(Grid->start, date_str, 5);
                        Grid->start[5] = '\0';
                    }
                    
                    Grid = Grid->next;
                }

                Grid = initial;
                while (Grid)
                {
                    for (i = 0; i <= Meteo->Seasons; i++)
                    {
                        Grid->twso[i] = 0.0;
                        Grid->length[i] = 0.0;
                        Grid->applied_n[i] = 0.0f;
                        Grid->cold_days[i] = 0;
                        Grid->crp->Seasons = 1;
                    }
                    Grid = Grid->next;
                }

                for (Day = 0; Day < Meteo->ntime; Day++)
                {
                    /* --- Check for new year to update fertilizer for this grid cell --- */
                    if (strlen(n_fert_file) > 0 && MeteoDay[Day] == 1) {
                        int current_year = MeteoYear[Day];
                        int year_index = current_year - Meteo->n_fert_start_year;
                    
                        if (year_index >= 0 && year_index < Meteo->n_fert_time_len) {
                            float annual_n_amount = Meteo->n_fertilizer_grid[year_index][Lat][Lon];
                    
                            int dekad = (int)Meteo->sowing_date_grid[Lat][Lon];
                            int sow_month = 1, sow_day = 1;
                            if (dekad >= 1 && dekad <= 36) {
                                sow_month = ((dekad - 1) / 3) + 1;
                                int subdek = ((dekad - 1) % 3) + 1;
                                sow_day = (subdek == 1) ? 1 : (subdek == 2) ? 11 : 21;
                            }
                    
                            SimUnit *tempGrid = initial;
                            while (tempGrid) {
                                if (tempGrid->mng->N_Fert_table != NULL) {
                                    tempGrid->mng->N_Fert_table->month = sow_month;
                                    tempGrid->mng->N_Fert_table->day = sow_day;
                                    tempGrid->mng->N_Fert_table->amount = annual_n_amount;
                    
                                    TABLE_D *next_app = tempGrid->mng->N_Fert_table->next;
                                    while (next_app != NULL) {
                                        next_app->amount = 0.0;
                                        next_app = next_app->next;
                                    }

                                    int season_idx = (current_year - Meteo->StartYear) + 1;
                                    if (season_idx >= 1 && season_idx <= Meteo->Seasons) {
                                        SimUnit *g = initial;
                                        while (g) {
                                            g->applied_n[season_idx] = annual_n_amount;
                                            g = g->next;
                                        }
                                    }
                                }
                                tempGrid = tempGrid->next;
                            }
                        }
                    }

                    Grid = initial;

                    /* Set the date struct */ 
                    memset(&current_date, 0, sizeof(current_date));
                    current_date.tm_year = MeteoYear[Day] - 1900;
                    current_date.tm_mday = 0 + MeteoDay[Day];
                    mktime(&current_date);

                    while (Grid)
                    {
                        /* Get data, states and rates from the Grid structure */
                        Crop = Grid->crp;
                        WatBal = Grid->soil;
                        Mng = Grid->mng;
                        Site = Grid->ste;

                        Emergence = Grid->emergence;

                        Temp = 0.5 * (Tmax[Lon][Lat][Day] + Tmin[Lon][Lat][Day]);
                        DayTemp = 0.5 * (Tmax[Lon][Lat][Day] + Temp);

                        /* Only simulate between start and end year */ 
                        if ((MeteoYear[Day] >= Meteo->StartYear && MeteoYear[Day] <= Meteo->EndYear) && (Meteo->Seasons >= Crop->Seasons))
                        {
                            IfSowing(Grid->start);

                            if (Crop->Sowing >= 1 && Crop->Emergence == 0)
                            {
                                if (EmergenceCrop(Emergence))
                                {
                                    InitializeCrop();
                                    InitializeWatBal();
                                    InitializeNutrients();
                                }
                            }

                            if (Crop->Sowing >= 1 && Crop->Emergence == 1)
                            {
                                if (Crop->st.Development <= (Crop->prm.DevelopStageHarvest) && Crop->GrowthDay < CycleLength)
                                {
                                    Astro();
                                    CalcPenman();
                                    CalcPenmanMonteith();

                                    EvapTra();

                                    RatesToZero();

                                    RateCalulationWatBal();
                                    Partioning();
                                    RateCalcultionNutrients();
                                    RateCalculationCrop();

                                    Crop->st.LAI = LeaveAreaIndex();

                                    IntegrationCrop();

                                    if (Temp >= 0.0f && Temp <= 7.0f && Crop->st.Development <= 0.3f) {
                                        int year_idx = MeteoYear[Day] - Meteo->StartYear + 1;
                                        if (year_idx >= 1 && year_idx <= Meteo->Seasons) {
                                            Grid->cold_days[year_idx] += 1;
                                        }
                                    }

                                    IntegrationWatBal();
                                    IntegrationNutrients();

                                    Crop->GrowthDay++;
                                }
                                else
                                {
                                    Grid->twso[Crop->Seasons] = Crop->st.storage;
                                    Grid->n_storage[Crop->Seasons] = Crop->N_st.storage;    // Final N storage at harvest
                                    Grid->p_storage[Crop->Seasons] = Crop->P_st.storage;    // Final P storage at harvest  
                                    Grid->k_storage[Crop->Seasons] = Crop->K_rt.storage;    // Final K storage at harvest
                                    Grid->length[Crop->Seasons] = Crop->GrowthDay;

                                    if (Meteo->Seasons == Crop->Seasons)
                                    {
                                        Output(files[Grid->file]);
                                        WriteOutputToNetCDF(&nc_output);
                                    }

                                    while (Crop->LeaveProperties != NULL)
                                    {
                                        wipe = Crop->LeaveProperties;
                                        Crop->LeaveProperties = Crop->LeaveProperties->next;
                                        free(wipe); 
                                    }

                                    Emergence = 0;
                                    Crop->TSumEmergence = 0;
                                    Crop->Emergence = 0;
                                    Crop->Sowing = 0;
                                    Crop->Seasons++;
                                }
                            }
                        }

                        Grid->crp = Crop;
                        Grid->soil = WatBal;
                        Grid->mng = Mng;
                        Grid->ste = Site;
                        Grid = Grid->next;
                    }
                }
            }
        }

        head = Meteo;
        Meteo = Meteo->next;
        CleanGridData(head); 

        if (strlen(n_fert_file) > 0) {
            CleanFertilizerData(head);
        }
        CleanMeteo(head); 
        free(head);
    }

    /* Close NetCDF file */
    printf("Closing NetCDF file...\n");
    CloseNetCDF(&nc_output);

    /* Close the output files and free the allocated memory */
    fptr = NULL;
    Grid = initial;
    while (Grid)
    {
        if (files[Grid->file] != fptr)
        {
            if (Grid->file < NumberOfFiles)
            {
                fptr = files[Grid->file];
                fclose(files[Grid->file]); 
            }
        }
        Grid = Grid->next;
    }

    Grid = initial;
    Clean(Grid); 
    free(files);

usage_error:
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s <sim_list> <meteolist>                          # Default TSUM and sowing\n", argv[0]);
    fprintf(stderr, "  %s <sim_list> <meteolist> <grid_data> <t1> <t2> <sow> [flags]  # Full gridded TSUM\n", argv[0]);
    fprintf(stderr, "  %s <sim_list> <meteolist> <grid_data> <sow> [flags]             # Default TSUM + gridded sowing\n", argv[0]);
    fprintf(stderr, "  %s <sim_list> <meteolist> <sow> [flags]                        # Default TSUM + sowing var name\n", argv[0]);
    fprintf(stderr, "\nOptional flags:\n");
    fprintf(stderr, "  --use-potential-nutrients\n");
    fprintf(stderr, "  --use-potential-evtra\n");
    fprintf(stderr, "  --n-fertilizer-nc <file.nc>\n");
    fprintf(stderr, "  --tsum1-offset <float>\n");
    fprintf(stderr, "  --tsum2-offset <float>\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  %s list.txt meteolist.txt                          # All defaults\n", argv[0]);
    fprintf(stderr, "  %s list.txt meteolist.txt all_griddata.nc tsum1 tsum2 sow_e1  # Full gridded\n", argv[0]);
    fprintf(stderr, "  %s list.txt meteolist.txt all_griddata.nc sow_e1 --tsum1-offset 50  # Default TSUM + gridded sowing\n", argv[0]);
    fprintf(stderr, "  %s list.txt meteolist.txt sow_e1 --tsum1-offset 50              # Default TSUM + sowing var name\n", argv[0]);
    exit(1);

    return 1;
}