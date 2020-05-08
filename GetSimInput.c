#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wofost.h"
#include "extern.h"

/* ----------------------------------------------------------------------------*/
/*  function GetMeteoInput()                                                   */
/*  Purpose: Get the names of the crop, soil, management and site files as     */
/*           as well if the simulation starting day and whether the simulation */
/*           starts at emergence (1) or sowing (0)                             */
/*-----------------------------------------------------------------------------*/

int GetSimInput(char *list)
{
    FILE *ifp;
     
    SimUnit *initial = NULL;
       
    int Emergence;
    int count;
  
    char line[MAX_STRING];
    char path[MAX_STRING];
    char cropfile[MAX_STRING];
    char soilfile[MAX_STRING];
    char sitefile[MAX_STRING];
    char management[MAX_STRING];
    char output[MAX_STRING];
    char start[MAX_STRING];
    char cf[MAX_STRING], sf[MAX_STRING], mf[MAX_STRING], site[MAX_STRING];
  
    
    ifp = fopen(list, "r");

    if (ifp == NULL) 
    {
        fprintf(stderr, "Can't open input, %s\n", list);
        exit(1);
    }
    
    count = 0;
    while(fgets(line, MAX_STRING, ifp)) {
        if(line[0] == '*' || line[0] == ' ' || line[0] == '\n'){
            continue;
        }
        
        sscanf(line,"%s %s %s %s %s %s %d %s" ,
            path, cf, sf, mf, site, start, &Emergence, output);
        
        memset(cropfile,'\0',MAX_STRING);
        memset(sitefile,'\0',MAX_STRING);
        memset(soilfile,'\0',MAX_STRING);
        memset(management,'\0',MAX_STRING);
                
        strncpy(cropfile, path, strlen(path));
        strncat(cropfile, cf, strlen(cf));

        strncpy(soilfile, path, strlen(path));
        strncat(soilfile, sf, strlen(sf));

        strncpy(management, path, strlen(path));
        strncat(management, mf, strlen(mf));

        strncpy(sitefile, path, strlen(path));
        strncat(sitefile, site, strlen(site));
        
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

        if (strlen(sf) >= MAX_STRING) exit(0);
        if (strlen(output) >= MAX_STRING) exit(0);    
        if (strlen(start) >= MAX_STRING) exit(0);  
        
        memset(Grid->output,'\0',MAX_STRING);
        memset(Grid->start,'\0',MAX_STRING);
        
        strncpy(Grid->output,output,strlen(output)); // Name og output file
        strncpy(Grid->start,start,strlen(start)); // Starting string month day of the simulations 
        
        Grid->file  = count++;            // number of elements in Grid carousel
        Grid->emergence = Emergence;      // Start the simulations at emergence (1) or at sowing (0)                
        Grid->crp->Sowing = 0;
        Grid->crp->Emergence = 0;         // Crop emergence has not yet occurred
        Grid->next = NULL;
    }
    
    fclose(ifp);
    
    /* Set Grid back to initial address */
    Grid = initial;
    
    /* Return the number of meteo files */
    return count;
}   