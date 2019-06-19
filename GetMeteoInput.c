#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wofost.h"
#include "extern.h"

void GetMeteoInput()
{
    FILE *ifp;
    
    Weather *initial = NULL;
    
    int StartYear;
    int EndYear;
    float lat;
    float lon;
    char filename[100];
    
    ifp = fopen("meteolist.txt", "r");

    if (ifp == NULL) 
    {
        fprintf(stderr, "Can't open input meteolist.txt\n");
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
            Meteo->next = malloc(sizeof(SimUnit));
            Meteo = Meteo->next;  
        }    
        strcpy(Meteo->file, filename);
        Meteo->StartYear = StartYear;
        Meteo->EndYear = EndYear;
        Meteo->lat = lat;
        Meteo->lon = lon;
        Meteo->next = NULL;
    }
        
    Meteo = initial;
    fclose(ifp);
  
}