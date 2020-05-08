#include "wofost.h"
#include "extern.h"

/* ---------------------------------------------------------------------------*/
/*  function List()                                                           */
/*  Purpose: Get the value of a user provided input table                     */
/* ---------------------------------------------------------------------------*/
     
float List(TABLE_D *Table)
{
     while (Table) 
    { 
        if (current_date.tm_mon == Table->month -1 &&
            current_date.tm_mday == Table->day && 
            MeteoYear[Day] <= Meteo->EndYear)
        {
            return Table->amount;
        }
        Table = Table->next;
    }
    
    return 0.;
}
