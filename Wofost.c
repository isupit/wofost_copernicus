#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"
#include "output_netcdf.h"
#include "input_griddata.h"
#include "input_fertilizer.h"

// --- Global variables for optional flags ---
// Initialize them to their default state (0 = off/false)
int use_potential_nutrients = 0;
int use_potential_evtra = 0;

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

    // We need at least 7 (program name + 6 mandatory args).
    // Allow up to 11 now.
    if (argc < 7 || argc > 11) {
        fprintf(stderr, "Usage: %s <sim_list> <meteo_list> <grid_data> <tsum1_var> <tsum2_var> <sow_var> [--use-potential-nutrients] [--use-potential-evtra] [--n-fertilizer-nc <file.nc>]\n", argv[0]);
        fprintf(stderr, "Example: %s list.txt meteolist.txt all_griddata.nc avg_tsum1_e1e1 avg_tsum2_e1e1 sow_e1 --n-fertilizer-nc fertilizer.nc\n", argv[0]);
        exit(0);
    }

    // These checks can remain the same
    if (strlen(argv[1]) >= MAX_STRING) exit(0);
    if (strlen(argv[2]) >= MAX_STRING) exit(0);
    if (strlen(argv[3]) >= MAX_STRING) exit(0);
    if (strlen(argv[4]) >= MAX_STRING) exit(0);
    if (strlen(argv[5]) >= MAX_STRING) exit(0);
    if (strlen(argv[6]) >= MAX_STRING) exit(0);

    memset(list, '\0', MAX_STRING);
    memset(meteolist, '\0', MAX_STRING);
    memset(grid_data_file, '\0', MAX_STRING); 
    memset(tsum1_var, '\0', MAX_STRING); 
    memset(tsum2_var, '\0', MAX_STRING); 
    memset(sow_var, '\0', MAX_STRING);
    memset(n_fert_file, '\0', MAX_STRING);


    // --- Parse the 6 mandatory arguments first ---
    strncpy(list, argv[1], strlen(argv[1]));
    strncpy(meteolist, argv[2], strlen(argv[2]));
    strncpy(grid_data_file, argv[3], strlen(argv[3])); 
    strncpy(tsum1_var, argv[4], strlen(argv[4])); 
    strncpy(tsum2_var, argv[5], strlen(argv[5])); 
    strncpy(sow_var, argv[6], strlen(argv[6]));

    // --- Loop through the OPTIONAL arguments and check for your specific flags ---
    for (int i = 7; i < argc; i++) {
        if (strcmp(argv[i], "--use-potential-nutrients") == 0) {
            use_potential_nutrients = 1; // Set flag to true
        } else if (strcmp(argv[i], "--use-potential-evtra") == 0) {
            use_potential_evtra = 1; // Set flag to true
        } else if (strcmp(argv[i], "--n-fertilizer-nc") == 0) {
            if (i + 1 < argc) { // Make sure a filename is provided
                strncpy(n_fert_file, argv[i + 1], MAX_STRING - 1);
                i++; // Increment i to skip the filename in the next iteration
            } else {
                fprintf(stderr, "Error: --n-fertilizer-nc flag requires a filename.\n");
                exit(1);
            }
        } else {
            // If the argument is unknown, print an error and exit
            fprintf(stderr, "Error: Unknown optional argument '%s'\n", argv[i]);
            fprintf(stderr, "Usage: %s <sim_list> <meteo_list> <grid_data> <tsum1_var> <tsum2_var> <sow_var> [--use-potential-nutrients] [--use-potential-evtra] [--n-fertilizer-nc <file.nc>]\n", argv[0]);
            exit(1);
        }
    }

    // --- MODIFIED: You can now check the flags later in your code ---
    printf("\n--- Configuration Summary ---\n");
    printf("Mandatory arguments loaded successfully.\n");
    printf("Optional flag --use-potential-nutrients set: %s\n", use_potential_nutrients ? "Yes" : "No");
    printf("Optional flag --use-potential-evtra set: %s\n", use_potential_evtra ? "Yes" : "No");
    printf("Optional N Fertilizer NetCDF provided: %s\n", strlen(n_fert_file) > 0 ? n_fert_file : "No");
    printf("---------------------------\n\n");

    /* --- Construct dynamic output filename --- */
    /* e.g. for the pair tsum1_e1a1, tsum2_e1a1 and sow_e1, the output filename becomes: wofost_results_e1a1_sow_e1.nc */
    char *tsum_suffix = strrchr(tsum1_var, '_');
    if (tsum_suffix && strlen(tsum_suffix) > 1) {
        tsum_suffix++; // Skip the '_'
    } else {
        tsum_suffix = "default";
        fprintf(stderr, "Warning: Could not extract tsum suffix from %s, using default\n", tsum1_var);
    }
    if (strlen(sow_var) < 1) {
        fprintf(stderr, "Error: sow_var is empty\n");
        exit(1);
    }
    snprintf(output_file, MAX_STRING, "wofost_results_%s_%s.nc", tsum_suffix, sow_var);
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
        // Grid->twso = (float*) malloc((Meteo->EndYear - Meteo->StartYear + 1) * sizeof(float));
        for (i = 0; i <= Meteo->Seasons; i++)
        {
            Grid->twso[i] = 0.0;
            Grid->length[i] = 0.0;
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

        /* Load crop grid data after meteo/weather dimensions are known */
        GetGridData(Meteo, grid_data_file, tsum1_var, tsum2_var, sow_var);

        /* ---  Load fertilizer data if provided --- */
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
                    /* Update Tsum based on grid position */
                    Grid->crp->prm.TempSum1 = Meteo->tsum1_grid[Lat][Lon];
                    Grid->crp->prm.TempSum2 = Meteo->tsum2_grid[Lat][Lon];
                    
                    /* Convert dekad (float) to "MM-DD" string  */
                    int dekad = (int)Meteo->sowing_date_grid[Lat][Lon];
                    if (dekad < 1 || dekad > 36) {
                        // Default or error handling; using Jan 1 as fallback
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
                        Grid->crp->Seasons = 1;
                    }
                    Grid = Grid->next;
                }

                for (Day = 0; Day < Meteo->ntime; Day++)
                // assume that the series start January first 
                {

                    /* --- Check for new year to update fertilizer for this grid cell --- */
                    if (strlen(n_fert_file) > 0 && MeteoDay[Day] == 1) {
                        int current_year = MeteoYear[Day];
                        int year_index = current_year - Meteo->n_fert_start_year;
                    
                        if (year_index >= 0 && year_index < Meteo->n_fert_time_len) {
                            // Get the annual fertilizer amount for the current year and grid cell
                            float annual_n_amount = Meteo->n_fertilizer_grid[year_index][Lat][Lon];
                    
                            // Get the sowing date (month and day) for the current grid cell
                            int dekad = (int)Meteo->sowing_date_grid[Lat][Lon];
                            int sow_month = 1, sow_day = 1; // Default fallback
                            if (dekad >= 1 && dekad <= 36) {
                                sow_month = ((dekad - 1) / 3) + 1;
                                int subdek = ((dekad - 1) % 3) + 1;
                                sow_day = (subdek == 1) ? 1 : (subdek == 2) ? 11 : 21;
                            }
                    
                            // Loop through all SimUnits and update their management data
                            SimUnit *tempGrid = initial;
                            while (tempGrid) {
                                if (tempGrid->mng->N_Fert_table != NULL) {
                                    // --- Reprogram the FIRST application event ---
                                    // Set the date to the sowing date
                                    tempGrid->mng->N_Fert_table->month = sow_month;
                                    tempGrid->mng->N_Fert_table->day = sow_day;
                                    // Set the amount to the FULL annual total
                                    tempGrid->mng->N_Fert_table->amount = annual_n_amount;
                    
                                    // --- Disable all SUBSEQUENT application events ---
                                    TABLE_D *next_app = tempGrid->mng->N_Fert_table->next;
                                    while (next_app != NULL) {
                                        next_app->amount = 0.0;
                                        next_app = next_app->next;
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
                        /* Get data, states and rates from the Grid structure and */
                        /* put them in the place holders */
                        Crop = Grid->crp;
                        WatBal = Grid->soil;
                        Mng = Grid->mng;
                        Site = Grid->ste;

                        // Start     = Grid->start;
                        Emergence = Grid->emergence; /* Start simulation at sowing or emergence */

                        Temp = 0.5 * (Tmax[Lon][Lat][Day] + Tmin[Lon][Lat][Day]);
                        DayTemp = 0.5 * (Tmax[Lon][Lat][Day] + Temp);

                        /* Only simulate between start and end year */ 
                        if ((MeteoYear[Day] >= Meteo->StartYear && MeteoYear[Day] <= Meteo->EndYear) && (Meteo->Seasons >= Crop->Seasons))
                        {
                            /* Determine if the sowing already has occurred */ 
                            IfSowing(Grid->start);

                            /* If sowing has occurred than determine the emergence */
                            if (Crop->Sowing >= 1 && Crop->Emergence == 0)
                            {
                                if (EmergenceCrop(Emergence))
                                {
                                    /* Initialize: set state variables */ 
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

                                    /* Calculate the evapotranspiration */ 
                                    EvapTra();

                                    /* Set the rate variables to zero */ 
                                    RatesToZero();

                                    /* Rate calculations */ 
                                    RateCalulationWatBal();
                                    Partioning();
                                    RateCalcultionNutrients();
                                    RateCalculationCrop();

                                    /* Write to the output files */ 
                                    // Output(output[Grid->file]);

                                    /* Calculate LAI */ 
                                    Crop->st.LAI = LeaveAreaIndex();

                                    /* State calculations */ 
                                    IntegrationCrop();
                                    IntegrationWatBal();
                                    IntegrationNutrients();

                                    /* Update the number of days that the crop has grown*/ 
                                    Crop->GrowthDay++;
                                }
                                else
                                {
                                    /* Write to the output files */ 
                                    Grid->twso[Crop->Seasons] = Crop->st.storage;
                                    Grid->length[Crop->Seasons] = Crop->GrowthDay;
                                    if (Meteo->Seasons == Crop->Seasons)
                                    {
                                        /* Original output to text file */
                                        Output(files[Grid->file]);

                                        WriteOutputToNetCDF(&nc_output);
                                    }

                                    /* Clean the LeaveProperties */ 
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

                        /* Store the daily calculations in the Grid structure */ 
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
        /* --- Clean up fertilizer data --- */
        if (strlen(n_fert_file) > 0) {
            CleanFertilizerData(head);
        }
        CleanMeteo(head); 
        free(head);
    }

    /* <<< 3. CLOSE NETCDF FILE >>> */
    printf("Closing NetCDF file...\n");
    CloseNetCDF(&nc_output);

    /* Return to the beginning of the list */ 
    Grid = initial;

    /* Close the output files and free the allocated memory */
    fptr = NULL;
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

    /* Go back to the beginning of the list */ 
    Grid = initial;
    Clean(Grid); 

    return 1;
}
