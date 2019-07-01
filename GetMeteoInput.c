#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wofost.h"
#include "extern.h"

/* --------------------------------------------------------------------------*/
/*  function GetMeteoInput()                                                 */
/*  Purpose: Get the names of meteo files and the start and end year of the  */
/*           of the simulations. Read the la and lon of the boxes. Store the */
/*           information of each cell in the Meteo struct.                   */
/* --------------------------------------------------------------------------*/

void GetMeteoInput(char *meteolist)
{
    FILE *ifp;
    
    Weather *initial = NULL;
    
    int StartYear;
    int EndYear;
    float lat;
    float lon;
    char filename[MAX_STRING];
    
    ifp = fopen(meteolist, "r");

    if (ifp == NULL) 
    {
        fprintf(stderr, "Can't open %s\n", meteolist);
        exit(1);
    }
    
    while (fscanf(ifp,"%s %d %d %f %f" , filename, &StartYear, &EndYear, &lat, &lon) != EOF) 
    {
        if (initial == NULL) 
        {
            Meteo = initial = malloc(sizeof(Weather));
        }
        else 
        {
            Meteo->next = malloc(sizeof(Weather));
            Meteo = Meteo->next;  
        }  
        
        if (strlen(filename) >= MAX_STRING) exit(0);
        
        memset(Meteo->file,'\0',MAX_STRING);
        strncpy(Meteo->file, filename, strlen(filename));
        
        Meteo->StartYear = StartYear;
        Meteo->EndYear = EndYear;
        Meteo->lat = lat;
        Meteo->lon = lon;
        Meteo->next = NULL;
    }
          
    Meteo = initial;
    fclose(ifp);
  
}