#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"
#include <time.h>

int main() {
    
    FILE **output;
    
    SimUnit *initial  = NULL;
    Weather *Meteo = NULL;
    Weather *head;
      
    int CycleLength   = 300;
    int Start;
    int Emergence;
    
    char name[100];
    
    /* Get the initial Grid address */
    initial = GetSimInput();    
    
    /* Get the meteo filenames and put them in the placeholder */
    GetMeteoInput();
    
    /* Allocate memory for the file pointers */
    output = malloc(sizeof(**output) * --count);
    
    /* Go back to the beginning of the list */
    Grid = initial;
    
    /* Open the output files */
    while (Grid)
    {   /* Make valgrind happy  */
        memset(name,0,100);
        
        memcpy(name, Grid->name, strlen(Grid->name)-4);
        snprintf(name, sizeof name, "%s%s%d%s", Grid->name, "-", Grid->file,".txt");
        
        output[Grid->file] = fopen(name, "w");
        header(output[Grid->file]);
        Grid = Grid->next;
    }
    
    /* Go back to the beginning of the list */
    Grid = initial;
    
    /* Set the starting year to all grid members */
    Grid->season = Meteo->StartYear;
    Grid = Grid->next;
    
    while (Meteo)
    {
        /* Get the meteodata */
        GetMeteoData(path, dateString, Meteo->Name);
    
        for (Day = 1; Day < 34333; Day++)
        {        
            /* Go back to the beginning of the list */
            Grid = initial;

            while (Grid)
            {
                /* Get data, states and rates from the Grid structure and */
                /* put them in the place holders */
                Crop      = Grid->crp;
                WatBal    = Grid->soil;
                Mng       = Grid->mng;
                Site      = Grid->ste;
                Start     = Grid->start;
                Emergence = Grid->emergence;
                
                Temp = 0.5 * (Tmax[Day] + Tmin[Day]);
                DayTemp = 0.5 * (Tmax[Day] + Temp);

                if (MeteoDay[Day] >= Start &&  MeteoYear[Day] == Grid->season)
                {
                    Crop->Sowing = 0;
                
                    /* Continue if the number of seasons has not been reached */
                    if ((Grid->season - Meteo->StartYear) <= Meteo->NumberOfYears)
                    {
                        Crop->Sowing = 1;
                    }
                }
 
                if (EmergenceCrop(Emergence))
                {                 
                    /* Initialize: set state variables */
                    InitializeCrop();
                    InitializeWatBal();
                    InitializeNutrients(); 
                }
                

                if (MeteoDay[Day] >= Start && Crop->Emergence == 1)
                {   
                    if (Crop->st.Development <= (Crop->prm.DevelopStageHarvest+0.05) && Crop->GrowthDay < CycleLength) 
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
                        Output(output[Grid->file]);   

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
                        Grid->season++;
                        Emergence = 0;
                        Crop->Emergence = 0;
                        Crop->Sowing    = 0;
                    }
                }

                /* Store the daily calculations in the Grid structure */
                Grid->crp  = Crop;
                Grid->soil = WatBal;
                Grid->mng  = Mng;
                Grid->ste  = Site;
                Grid = Grid->next;
            }

            /* Update time */
            simTime.tm_hour = 18;
            simTime.tm_mday++;
            mktime(&simTime);
        }   
    head = Meteo;
    Meteo = Meteo->next;
    free(head);
    }
    free(Meteo);
    
    /* Return to the beginning of the list */
    Grid = initial;
    
    /* Close the output files and free the allocated memory */
    while(Grid)
    {
        fclose(output[Grid->file]);
        Grid = Grid->next;
    }
    free(output);

    /* Go back to the beginning of the list */
    Grid = initial;
    Clean(Grid);

    return 1;
}
