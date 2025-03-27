#include <math.h>
#include "wofost.h"
#include "extern.h"

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

float ext_min(float x1, float x2, float x3)
{

    return(fmin(fmin(x1,x2),x3));
}

