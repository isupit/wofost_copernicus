#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wofost.h"
#include "extern.h"

void GetSimInput()
{
    FILE *ifp;
     
    SimUnit *initial = NULL;
       
    int Emergence;
    int Start;
    int count;
  
    char path[100];
    char cropfile[100];
    char soilfile[100];
    char sitefile[100];
    char management[100];
    char cf[100], sf[100], mf[100], site[100];
  
    
    ifp = fopen("list.txt", "r");

    if (ifp == NULL) 
    {
        fprintf(stderr, "Can't open input list.txt\n");
        exit(1);
    }
    
    count = 0;
    while (fscanf(ifp,"%7s %11s %7s %12s %12s %d %d" ,
            path, cf, sf, mf, site, &Start, &Emergence)
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
        Grid->file  = ++count;          // number of elements in Grid carousel
        strcpy(Grid->name,sf);          // Set the soil filename as ouput file name
        Grid->emergence = Emergence;    // Start the simulations at emergence (1) or at sowing (0)
        Grid->start = Start;             // Starting day of the simulations     
        Grid->crp->Sowing = 0;
        Grid->crp->Emergence = 0;       // Crop emergence has not yet occurred
        
        Grid->next = NULL;

    }
    
    fclose(ifp);
    
    /* Set Grid back to initial address */
    Grid = initial;
}   