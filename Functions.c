#include <stdio.h>
#include <ctype.h>
#include "wofost.h"
#include "extern.h"

/* ---------------------------------------------------------------------------*/
/*  function List()                                                           */
/*  Purpose: Get the value of a user provided input table                     */
/* ---------------------------------------------------------------------------*/

//float List_D(TABLE *Table, float *X)
//{
//if (*X < Table->x)  return 0.;
//
//while (Table->next) 
//{
//    if (*X == Table->x)  
//	return (Table->y);
//    Table = Table->next;
//}        
//
//return 0.;
//}       
  
float List(TABLE_D *Table)
{
    struct tm application_date = { 0 };
    struct tm current_date = { 0 };
    
    current_date.tm_year = MeteoYear[Day] -1900;
    current_date.tm_mday =  0 + MeteoDay[Day];
    mktime(&current_date);

    while (Table->next) 
    {
        application_date.tm_year = MeteoYear[Day] -1900;
        application_date.tm_mon  = Table->month -1;
        application_date.tm_mday = Table->day;
        mktime(&application_date);   
  
        if (application_date.tm_mon == current_date.tm_mon &&
            application_date.tm_mday== current_date.tm_mday && 
            MeteoYear[Day] <= Meteo->EndYear)
        {
            return Table->amount;
        }
        Table = Table->next;
    }
    
    return 0.;
}



float limit(float a, float b, float c)
{
    if (c < a) return a;
    else if (c >= a && c <= b)  return c;
    else return b;
   }


float notnul(float x)
{
    if (x != 0.) return x;
    else return 1.;
   }

float insw(float x1, float x2, float x3)
{
    return ((x1 < 0) ? x2 : x3);
}


int leap_year(int year)
{
    if ((year % 400 == 0) || ( ( year % 100 != 0) && (year % 4 == 0 )))
        return 366;
    else
        return 365;
}


float min(float a, float b)
{
    return ((a < b) ? a : b);
}


float max(float a, float b)
{
    return ((a > b) ? a : b);
}

