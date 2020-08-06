#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "astro.h"
#include "extern.h"
#include "wofost.h"

/* ---------------------------------------------------------------------*/
/*  function Astro()                                                    */
/*  Purpose: Calculation of the astronomical parameters used in Wofost  */
/*                                                                      */
/*  Originally written in Fortran by:                                   */
/*         Daniel van Kraalingen, April 1991                            */
/*         revised Allard de Wit, January 2011                          */
/* ---------------------------------------------------------------------*/

int Astro()
{
    float Declination;
    float AOB;
    
    float day_fl  = 1.0 + (float) current_date.tm_yday;
    float Latitude = lats[lat];
    
    if (fabs(Latitude) > 90.) return 0;  

    /* We start at Day= 1, we do not use Day = 0 */
    Declination    = -asin(sin(23.45*RAD)*cos(2.*PI*(day_fl+10.)/365.));
    SolarConstant  = 1367.*(1.+0.033*cos(2.*PI*(day_fl -10)/365.));
  
    SinLD = sin(RAD*Latitude)*sin(Declination);
    CosLD = cos(RAD*Latitude)*cos(Declination);
    AOB   = SinLD/CosLD;
    
   /* Astronomical day length */
    Daylength = max(0,min(24.,12.0*(1.+2.*asin(AOB)/PI)));
    
    /* SunRise-Set */
    SunRise = 12. - 0.5*Daylength;
    SunSet  = 12. + 0.5*Daylength;
    
    /* Photoactive day length */
    PARDaylength = max(0,min(24.,12.0*(1.+2.*asin((-sin(ANGLE*RAD)+SinLD)/CosLD)/PI)));
    
    /* Integrals of sine of solar height */
    if (AOB <= 1.0)  
        DSinBE = 3600.*(Daylength*(SinLD+0.4*(SinLD*SinLD + CosLD*CosLD*0.5))+
            12.*CosLD*(2.+3.*0.4*SinLD)*sqrt(1.-AOB*AOB)/PI);
    else
        DSinBE = 3600.*(Daylength*(SinLD+0.4*(SinLD*SinLD + CosLD*CosLD*0.5)));
   
    return 1;
}
