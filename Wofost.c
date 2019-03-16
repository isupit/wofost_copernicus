#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"
#include <time.h>

int main() {
    FILE *ifp;
    FILE **output;
    
    SimUnit *Grid = NULL;
    SimUnit *initial = NULL;
       
    int Emergence;
    int Start;
    int CycleLength   = 300;
    int count;
    
    char path[100];
    char cropfile[100];
    char soilfile[100];
    char sitefile[100];
    char management[100];
    char dateString [100];
    char place[15];
    char name[100];
      
    char cf[100], sf[100], mf[100], site[100];
  
    Step = 1.;    
    
    ifp = fopen("list.txt", "r");

    if (ifp == NULL) 
    {
        fprintf(stderr, "Can't open input list.txt\n");
        exit(1);
    }
    
    count = 0;
    while (fscanf(ifp,"%7s %11s %7s %12s %10s %10s %4s %d %d" ,
            path, cf, sf, mf, site, dateString, place, &Start, &Emergence)
            != EOF) 
    {    
        strncpy(cropfile, path, 98);
        strncat(cropfile, cf, 98);

        strncpy(soilfile, path, 98);
        strncat(soilfile, sf, 98);

        strncpy(management, path, 98);
        strncat(management, mf, 98);

        strncpy(sitefile, path, 98);
        strncat(sitefile, site, 98);
        
        /* count the number of output files */
        /* number is the index number of the list of file pointers */
        if (initial == NULL) 
        {
            Grid = initial =  malloc(sizeof(SimUnit));
        }
        else 
        {
            Grid->next = malloc(sizeof(SimUnit));
            Grid = Grid->next;  
        }
        
        GetCropData(Grid->crp   = malloc(sizeof(Plant)), cropfile); 
        GetSiteData(Grid->ste   = malloc(sizeof(Field)), sitefile);
        GetManagement(Grid->mng = malloc(sizeof(Management)), management);
        GetSoilData(Grid->soil  = malloc(sizeof(Soil)), soilfile);

        Grid->start = Start;            // Start day (=day number)
        Grid->file  = count++;          // number of elements in Grid carousel
        strcpy(Grid->name,cf);          // Crop file name
        Grid->emergence = Emergence;    // Start the simulations at emergence (1) or at sowing (0)
        Grid->next = NULL;

    }
    
    /* Close the input file */
    fclose(ifp);
    
    /* Set Grid back to the initial address */
    Grid = initial;   
    
    /* Allocate memory for the file pointers */
    output = malloc(sizeof(**output) * --count);
    
    /* Open the output files */
    while (Grid)
    {   /* Make valgrind happy  */
        memset(name,0,100);
        
        memcpy(name, Grid->name, strlen(Grid->name)-4);
        snprintf(name, sizeof name, "%s%s%d%s", Grid->name, "-", Grid->file,".txt");
        printf("%20s\n",Grid->name);
        output[Grid->file] = fopen(name, "w");
        header(output[Grid->file]);
        Grid = Grid->next;
    }
    
    /* Get the meteodata */
    GetMeteoData(path, dateString, place);
    
    /* Set Grid back to the initial address */
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
        printf("%20s\n",Grid->name);
        
        for (Day = 1; Day < 11348; Day++)
        {        
            /* Go back to the beginning of the list */
            Grid = initial;

            Temp = 0.5 * (Tmax[Day] + Tmin[Day]);
            DayTemp = 0.5 * (Tmax[Day] + Temp);

            Astro();
            CalcPenman();
            
            if (Day >= Start && Crop->Emergence == 0)
            {
                if (EmergenceCrop(Emergence))
                {                 
                    /* Initialize: set state variables */
                    InitializeCrop();
                    InitializeWatBal();
                    InitializeNutrients(); 
                }
            }
            
            if (Day >= Start && Crop->Emergence == 1)
            {   
                if (Crop->st.Development <= Crop->prm.DevelopStageHarvest && Crop->GrowthDay < CycleLength) 
                {
                   /* Calculate the evapotranspiration */
                    EvapTra();
                    
                    /* Set the rate variables to zero */
                    RatesToZero();
                    
                     /* Rate calculations */
                    RateCalulationWatBal();
                    RateCalcultionNutrients();
                    RateCalculationCrop();
                    
                    /* Calculate LAI */
                    Crop->st.LAI = LeaveAreaIndex();         
                                        
                    /* State calculations */
                    IntegrationWatBal();
                    IntegrationNutrients();
                    IntegrationCrop();
                    
                    /* Update the number of days that the crop has grown*/
                    Crop->GrowthDay++;
                }
                /* Crop development finished or harvest occurred */
                else
                {
                    /* Write to the output files */
                    Output(output[Grid->file]); 
                    
                    /* Start the next year */
                    Crop->Emergence = 0;
                    Start += leap_year(Year[Day]);        
                }
            }
            /* Update time */
            simTime.tm_mday++;
            mktime(&simTime);
          
        } 
        
        /* Close the output files */
        fclose(output[Grid->file]);
           
        /* Store the daily calculations in the Grid structure */
        Grid->crp  = Crop;
        Grid->soil = WatBal;
        Grid->mng  = Mng;
        Grid->ste  = Site;
        Grid = Grid->next;
    }    
    
    /* Return to the beginning of the list */
    Grid = initial;
    
    //free(output);

    /* Go back to the beginning of the list */
    Grid = initial;
    //Clean(Grid);

    return 1;
}


