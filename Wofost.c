#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"
#include "output_netcdf.h"
#include "input_griddata.h"   /* <<< ADDED: Include for TSM input functions */

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
    char grid_data_file[MAX_STRING]; // +++ ADD THIS LINE TO DECLARE THE VARIABLE
    char tsum1_var[MAX_STRING]; // +++ ADDED: For dynamic tsum1 variable name
    char tsum2_var[MAX_STRING]; // +++ ADDED: For dynamic tsum2 variable name
    char sow_var[MAX_STRING];   // +++ ADDED: For dynamic sowing variable name
    char output_file[MAX_STRING]; // +++ ADDED: For dynamic output filename

    Step = 1.; 

    // MODIFIED: Check for 7 arguments now
    if (argc != 7) {
        fprintf(stderr, "Usage: %s <sim_list> <meteo_list> <grid_data_netcdf_file> <tsum1_var> <tsum2_var> <sow_var>\n", argv[0]);
        fprintf(stderr, "Example: %s list.txt meteolist.txt all_griddata_cropped.nc avg_tsum1_e1e1 avg_tsum2_e1e1 sow_e1\n", argv[0]);
        exit(0);
    }
    if (strlen(argv[1]) >= MAX_STRING) exit(0);
    if (strlen(argv[2]) >= MAX_STRING) exit(0);
    if (strlen(argv[3]) >= MAX_STRING) exit(0); // <<< ADDED: Check for new argument
    if (strlen(argv[4]) >= MAX_STRING) exit(0); // <<< ADDED: Check for tsum1_var
    if (strlen(argv[5]) >= MAX_STRING) exit(0); // <<< ADDED: Check for tsum2_var
    if (strlen(argv[6]) >= MAX_STRING) exit(0); // <<< ADDED: Check for sow_var

    memset(list, '\0', MAX_STRING);
    memset(meteolist, '\0', MAX_STRING); // empty the memory string
    memset(grid_data_file, '\0', MAX_STRING); // <<< ADDED: Clear tsum_file string
    memset(tsum1_var, '\0', MAX_STRING); // <<< ADDED
    memset(tsum2_var, '\0', MAX_STRING); // <<< ADDED
    memset(sow_var, '\0', MAX_STRING);   // <<< ADDED

    strncpy(list, argv[1], strlen(argv[1]));
    strncpy(meteolist, argv[2], strlen(argv[2]));
    strncpy(grid_data_file, argv[3], strlen(argv[3])); // +++ ADDED: Copy the filename
    strncpy(tsum1_var, argv[4], strlen(argv[4])); // <<< ADDED
    strncpy(tsum2_var, argv[5], strlen(argv[5])); // <<< ADDED
    strncpy(sow_var, argv[6], strlen(argv[6]));   // <<< ADDED

    /* Construct dynamic output filename */
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

    /* <<< 1. SETUP NETCDF FILE >>> */
    NcFile nc_output;
    printf("Setting up NetCDF output file '%s'...\n", output_file);
    SetupNetCDF(output_file, &nc_output, Meteo->nlat, Meteo->nlon, tsum1_var, tsum2_var, sow_var);

    while (Meteo)
    {
        /* Get the meteodata */ 
        if (GetMeteoData(Meteo) != 1)
        {
            fprintf(stderr, "Cannot get meteo data.\n");
            exit(0);
        }

        /* <<< MODIFIED: Load ALL grid data AFTER meteo dimensions are known >>> */
        GetGridData(Meteo, grid_data_file, tsum1_var, tsum2_var, sow_var);

        printf("running %d - %d\n", Meteo->StartYear, Meteo->EndYear);

        for (Lon = 0; Lon < Meteo->nlon; Lon++)
        {
            for (Lat = 0; Lat < Meteo->nlat; Lat++)
            {
                if (Mask[Lon][Lat] != 1) 
                {
                    continue;
                }

                /* <<< ADDED: Update TSM values for this grid cell >>> */
                Grid = initial;
                while(Grid) {
                    // This part updates Tsum (existing logic)
                    Grid->crp->prm.TempSum1 = Meteo->tsum1_grid[Lat][Lon];
                    Grid->crp->prm.TempSum2 = Meteo->tsum2_grid[Lat][Lon];
                    
                    // +++ MODIFIED: Convert dekad (float) to "MM-DD" string and set emergence flag
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
                    Grid->emergence = 1; // Hardcoded based on original list.txt example
                
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
                                        /* ORIGINAL OUTPUT TO TEXT FILE */
                                        Output(files[Grid->file]);

                                        /* <<< 2. WRITE TO NETCDF FILE (SIMPLIFIED CALL) >>> */
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
        CleanGridData(head); // <<< MODIFIED: Replaces CleanTSumGrids
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
