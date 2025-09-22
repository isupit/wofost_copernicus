#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
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

// FIXED: Initialize all variables to safe defaults
    use_gridded_tsum = 0;
    memset(list, '\0', MAX_STRING);
    memset(meteolist, '\0', MAX_STRING);
    memset(grid_data_file, '\0', MAX_STRING); 
    memset(tsum1_var, '\0', MAX_STRING); 
    memset(tsum2_var, '\0', MAX_STRING); 
    memset(sow_var, '\0', MAX_STRING);
    memset(n_fert_file, '\0', MAX_STRING);

    // FIXED: Improved argument parsing that handles flags in BOTH modes
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

int arg_index = 3;  // Start parsing from argv[3]

    // FIXED: Check if we have gridded TSUM mode (exactly 6 positional args before flags)
    use_gridded_tsum = 0;
    if (argc >= 7) {  // 2 mandatory + 4 gridded = 6 positional args minimum
        // Check if argv[3] looks like a filename (not a flag)
        if (strncmp(argv[3], "--", 2) != 0) {
            // Check we have all 4 gridded arguments
            if (argc < 7) {
                fprintf(stderr, "Error: Gridded mode requires 6 positional arguments\n");
                goto usage_error;
            }
            
            // Check string lengths for gridded arguments
            if (strlen(argv[3]) >= MAX_STRING || strlen(argv[4]) >= MAX_STRING ||
                strlen(argv[5]) >= MAX_STRING || strlen(argv[6]) >= MAX_STRING) {
                fprintf(stderr, "Error: One of the gridded arguments is too long\n");
                exit(1);
            }
            
            // Parse the 4 gridded arguments
            strncpy(grid_data_file, argv[3], strlen(argv[3])); 
            strncpy(tsum1_var, argv[4], strlen(argv[4])); 
            strncpy(tsum2_var, argv[5], strlen(argv[5])); 
            strncpy(sow_var, argv[6], strlen(argv[6]));
            use_gridded_tsum = 1;
            
            arg_index = 7;  // Next args are optional flags
        } else {
            // argv[3] starts with "--", so it's a flag, use default mode
            use_gridded_tsum = 0;
            arg_index = 3;
        }
    } else {
        // Less than 7 arguments, definitely default mode
        use_gridded_tsum = 0;
        arg_index = 3;  // Start parsing flags from argv[3]
    }

    // FIXED: Parse optional flags from arg_index onwards (works for BOTH modes)
    for (int i = arg_index; i < argc; i++) {
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
            fprintf(stderr, "Error: Unknown optional argument '%s'\n", argv[i]);
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
    if (use_gridded_tsum) {
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
        } else {
            tsum_suffix = "default";
            fprintf(stderr, "Warning: Could not extract tsum suffix from %s, using default\n", tsum1_var);
        }
        
        // FIXED: Safer string construction
        char suffix_part[MAX_STRING/2];
        snprintf(suffix_part, sizeof(suffix_part), "%s_%s", tsum_suffix, sow_var);
        snprintf(output_file, MAX_STRING, "wofost_results_%s.nc", suffix_part);
    } else {
        // FIXED: Include offset info in default mode filename
        if (tsum1_offset != 0.0f || tsum2_offset != 0.0f) {
            snprintf(output_file, MAX_STRING, "wofost_results_default_t1+%.0f_t2+%.0f.nc", 
                     tsum1_offset, tsum2_offset);
        } else {
            snprintf(output_file, MAX_STRING, "wofost_results_default.nc");
        }
    }
    
    printf("Constructed output file: %s\n", output_file);

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
    SetupNetCDF(output_file, &nc_output, Meteo->nlat, Meteo->nlon, Meteo->Seasons, tsum1_var, tsum2_var, sow_var);

    while (Meteo)
    {
        /* Get the meteodata */ 
        if (GetMeteoData(Meteo) != 1)
        {
            fprintf(stderr, "Cannot get meteo data.\n");
            exit(0);
        }
    
        /* Load crop grid data conditionally */
        if (use_gridded_tsum) {
            printf("Using gridded TSUM1/TSUM2 data from %s\n", grid_data_file);
            GetGridData(Meteo, grid_data_file, tsum1_var, tsum2_var, sow_var);
        } else {
            printf("Using default TSUM1/TSUM2 values from crop parameter file\n");
            SetDefaultTSUM(Meteo);
            if (strlen(sow_var) > 0) {
                LoadSowingDateOnly(Meteo, grid_data_file, sow_var);
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
                    
                    /* Convert dekad (float) to "MM-DD" string  */
                    int dekad = (int)Meteo->sowing_date_grid[Lat][Lon];
                    if (dekad < 1 || dekad > 36) {
                        strncpy(Grid->start, "01-01", 5);
                        Grid->start[5] = '\0';
                    } else {
                        int month = ((dekad - 1) / 3) + 1;
                        int subdek = ((dekad - 1) % 3) + 1;
                        int day = (subdek == 1) ? 1 : (subdek == 2) ? 11 : 21;
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
    fprintf(stderr, "  %s <sim_list> <meteolist> [optional flags]                           # Default TSUM values\n", argv[0]);
    fprintf(stderr, "  %s <sim_list> <meteolist> <grid_data> <t1> <t2> <sow> [optional flags]  # Gridded TSUM\n", argv[0]);
    fprintf(stderr, "\nOptional flags:\n");
    fprintf(stderr, "  --use-potential-nutrients\n");
    fprintf(stderr, "  --use-potential-evtra\n");
    fprintf(stderr, "  --n-fertilizer-nc <file.nc>\n");
    fprintf(stderr, "  --tsum1-offset <float>\n");
    fprintf(stderr, "  --tsum2-offset <float>\n");
    fprintf(stderr, "\nExamples:\n");
    fprintf(stderr, "  %s list.txt meteolist.txt --tsum1-offset 50 --tsum2-offset -20          # Default TSUM + offsets\n", argv[0]);
    fprintf(stderr, "  %s list.txt meteolist.txt all_griddata.nc tsum1 tsum2 sow --tsum1-offset 50  # Gridded TSUM + offset\n", argv[0]);
    exit(1);



    return 1;
}