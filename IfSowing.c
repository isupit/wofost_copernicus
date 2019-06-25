#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

void IfSowing(char* dateString)
{
    struct tm sowing_date;
    struct tm current_date;
     
    int month, start_day;
     
    if(sscanf(dateString, "%d-%d", &month, &start_day) != EOF)
    {
        sowing_date.tm_year = MeteoYear[Day] -1900;
        sowing_date.tm_mon = month - 1;
        sowing_date.tm_mday = start_day;
        sowing_date.tm_hour = 0;
        sowing_date.tm_min = 0;
        sowing_date.tm_sec = 0;
    }
    
    current_date.tm_year = MeteoYear[Day] -1900;
    current_date.tm_yday = MeteoDay[Day];
    current_date.tm_hour = 0;
    current_date.tm_min = 0;
    current_date.tm_sec = 0;
    
    if (difftime(mktime(&sowing_date), mktime(&current_date)) == 0 && 
            MeteoYear[Day] <= Meteo->EndYear)
    {
        Crop->Sowing = 1;
    }
    
}