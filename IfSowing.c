#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "wofost.h"
#include "extern.h"

/* ---------------------------------------------------------------------------*/
/*  function void IfSowing    ()                                              */
/*  Purpose: Checks whether sowing has occurred. Note that if the emergence   */
/*           flag is set to 1 the crop simulation starts the next day. If it  */
/*           is set to 0 the Emergence date has to be established.            */
/* ---------------------------------------------------------------------------*/

void IfSowing(char* dateString)
{
    struct tm sowing_date = { 0 };
    struct tm current_date = { 0 };
    int month, start_day;
     
    if(sscanf(dateString, "%d-%d", &month, &start_day) != EOF)
    {
        sowing_date.tm_year = MeteoYear[Day] -1900;
        sowing_date.tm_mon  = month-1;
        sowing_date.tm_mday = start_day;
        mktime(&sowing_date);   
    }
    
    current_date.tm_year = MeteoYear[Day] -1900;
    current_date.tm_mday =  0 + MeteoDay[Day];
    mktime(&current_date);
    
    if (sowing_date.tm_mon == current_date.tm_mon &&
        sowing_date.tm_mday== current_date.tm_mday && 
        MeteoYear[Day] <= Meteo->EndYear)
    {
        Crop->Sowing = 1;
    }
}