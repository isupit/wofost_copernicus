/* original: Penman.for from I.G.A.M. Noy and C.A. van Diepen, 
   van Kraalingen, and Allard de Wit, Sep 2011
*/

#include <stdio.h>
#include <math.h>
#include "astro.h"
#include "penman.h"
#include "wofost.h"
#include "extern.h"
#include "assim.h"

/* ---------------------------------------------------------------------*/
/*  function CalcPenman()                                               */
/*  Purpose: Calculation of ETO evapotranspiration     mm d-1           */
/*                          ESO soil evaporation       mm d-1           */
/*                          EO  open water evaporation mm d-1           */
/*                                                                      */
/*     Originally written in Fortran by:                                */
/*         I.G.A.M. Noy and C.A. van Diepen, September 1986             */
/*         revised van Kraalingen, April, van Diepen, October 1991      */
/*         revised van Kraalingen, and Allard de Wit, Sep 2011          */
/* ---------------------------------------------------------------------*/

void CalcPenman()
{
    float RelSunShineDuration;
    //float Tmpa;
    float Tdif;
    float BU;
    float Pbar;
    float Gamma;
    float Ea;
    //float Eac;
    float delta;
    float RB;
    //float Rnc;
    float Rnw; 
    float Rns; 
    float VapourP; 
    float SaturatedVap;
    
    float Psycon = 0.67;    // psychrometric instrument constant (mbar/Celsius-1)
    float Refcfw = 0.05;    // albedo for a water surface                        
    float Refcfs = 0.15;    // albedo for a soil surface                         
    //float Refcfc = 0.25;    // albedo for a  canopy                              
    float Lhvap  = 2.45e6;  // latent heat of evaporation of water (J/kg=J/mm)  
    float Stbc   = 4.9e-3;  // Stefan Boltzmann constant (J/m2/d/K4) */
    /* Preparatory calculations: mean daily temperature, temperature difference */
    /* (Celsius) and the Bu coefficient Bu of the wind function (depends  on    */ 
    /* temperature difference)                                                  */
    
    //Tmpa  = (Tmin[Day] + Tmax[Day])/2.;
    Tdif  = Tmax[Day][lat][lon] - Tmin[Day][lat][lon];
    BU    = 0.54 + 0.35 * limit(0.,1.,(Tdif-12.)/4.);

    /* Barometric pressure (mbar)             */
    /* Psychrometric constant (mbar/Celsius)  */
    Pbar  = 1013.*exp (-0.034*Altitude/(Temp + 273.));
    Gamma = Psycon * Pbar/1013.;


    /* Saturated vapour pressure according to equation of Goudriaan     */
    /* (1977) derivative of SVAP with respect to temperature, i.e.      */
    /* slope of the SVAP-temperature curve (mbar/Celsius).              */
            
    /* Measured vapour pressure not to exceed saturated vapour pressure */

    SaturatedVap  = 6.10588 * exp(17.32491 * Temp/(Temp+238.102));
    delta         = 238.102 * 17.32491 * SaturatedVap/pow((Temp +238.102),2);
    VapourP       = min(Vapour[Day][lat][lon],SaturatedVap);

    /* The expression n/N (RelLSSD) from the Penman formula is estimated   */
    /* from the Angstrom formula: RI=RA(A+B.n/N) -> n/N=(RI/RA-A)/B,       */
    /* where RI/RA is the atmospheric transmission obtained by a CALL      */
    /* to ASTRO: */
              
    RelSunShineDuration = limit(0.,1.,(AtmosphTransm-AngstA)/AngstB);

    /* Terms in Penman formula, for water, soil and canopy            */
    /* Net outgoing long-wave radiation (J/m2/d) acc. to Brunt (1932) */
    RB  = Stbc * pow((Temp+273.),4) * (0.56-0.079 * sqrt(VapourP)) *
              (0.1 + 0.9 * RelSunShineDuration);

    /* Net absorbed radiation, expressed in mm/d */
    Rnw = (Radiation[Day][lat][lon] * (1.-Refcfw)-RB)/Lhvap;
    Rns = (Radiation[Day][lat][lon] * (1.-Refcfs)-RB)/Lhvap;
    //Rnc = (Radiation[Day] * (1.-Refcfc)-RB)/Lhvap;

    /* Evaporative demand of the atmosphere (mm/d)  */
    Ea  = 0.26 * max (0.,(SaturatedVap-VapourP)) * (0.5+BU * Windspeed[Day][lat][lon]);
    //Eac = 0.26 * max (0.,(SaturatedVap-VapourP)) * (1.0+BU * Windspeed[Day]);
   
    /* Penman formula (1948)                */
    /* Ensure reference evaporation >= 0.   */
    /* Convert to cm/day                    */
    Penman.E0  = max(0., 0.1 * (delta*Rnw + Gamma*Ea)/(delta + Gamma));
    Penman.ES0 = max(0., 0.1 * (delta*Rns + Gamma*Ea)/(delta + Gamma));
    //Penman.ET0 = max(0., 0.1 * (delta*Rnc + Gamma*Eac)/(delta + Gamma));
    
}

void CalcPenmanMonteith(float Frac, float R_turb, float RB_water, float RB_heat, float R_stomata, float Radiation, float *RnetAbs,  float *PT)
{
    float Temperature;
    float Rnl_Tmp;
    float CskyRad;
    float Psr, Ptr, Ptd;
    float BlackBRad, RnetOut;
    
    Temperature = DTemp + Dif;

    // Clear Sky radiation [J/m2/DAY] from Angot TOA radiation
    // the latter is found through a call to astro()
    CskyRad = (0.75 + (2e-05 * Altitude)) * AngotRadiation;
    
    if (CskyRad > 0)
    {
        BlackBRad  = BOLTZM * pow((Temperature +273.),4); // Black body radiation 
        Rnl_Tmp  = BlackBRad * (0.56 -0.079 * sqrt(10*Vapour[Day][lat][lon])); //Note that Vap here is hPa!
 
        // net absorbed radiation sunlit
        RnetOut = Rnl_Tmp * (0.1+0.9*CskyRad) * Frac; // Net outgoing radiation
        *RnetAbs = Radiation - RnetOut;  // Net absorbed sunlit radiation

        // Intermediate variable related to resistances
        Psr    = PSYCH * (RB_water + R_turb + R_stomata)/(RB_heat + R_turb);

        // Radiation-determined term
        Ptr    = *RnetAbs * Delta /(Delta + Psr)/LHVAP;

        // Vapour pressure-determined term
        Ptd    = (VHCA * VapDeficit/(RB_heat + R_turb))/(Delta + Psr)/LHVAP;

        // Potential evaporation or transpiration
        *PT     = max(1.e-10,Ptr + Ptd);
    }
    else
    {
        Penman.ET0 = 0.;      
    }    
}