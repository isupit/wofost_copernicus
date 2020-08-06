/* original: Penman.for from I.G.A.M. Noy and C.A. van Diepen, */
/* van Kraalingen, and Allard de Wit, Sep 2011                 */

#include <stdio.h>
#include <math.h>
#include "penman.h"
#include "wofost.h"
#include "extern.h"

/* -----------------------------------------------------------------*/
/*  function EvapTra()                                              */
/*  Purpose: Calculates the water stress and the transpiration rate */
/* -----------------------------------------------------------------*/     
void EvapTra() {   
    float KDiffuse;
    float MaxReductionOxygenStress;
    float OxygenStress;
    float SoilMoistureAeration;
          
    KDiffuse = Afgen(Crop->prm.KDiffuseTb, &(Crop->st.Development));      
    Evtra.MaxEvapWater = Penman.E0 * exp(-0.75 * KDiffuse * Crop->st.LAI);
    Evtra.MaxEvapSoil  = max(0., Penman.ES0 * exp(-0.75 * KDiffuse * Crop->st.LAI));
    
    if (Crop->prm.Airducts) 
    {
        /* Critical soil moisture content for aeration */
        SoilMoistureAeration = WatBal->ct.MoistureSAT - WatBal->ct.CriticalSoilAirC;
        
        /* Count days since start oxygen shortage (up to 4 days) */
        if (WatBal->st.Moisture >= SoilMoistureAeration) {
            Crop->DaysOxygenStress = min(Crop->DaysOxygenStress++, 4.);
        }
        else 
        {
            Crop->DaysOxygenStress = 0.;
        }
        
        /* Maximum reduction reached after 4 days */
        MaxReductionOxygenStress = limit (0.,1.,(WatBal->ct.MoistureSAT - WatBal->st.Moisture)/
                (WatBal->ct.MoistureSAT - SoilMoistureAeration));
        
        OxygenStress   = MaxReductionOxygenStress + 
                (1.-Crop->DaysOxygenStress/4.)*(1.-MaxReductionOxygenStress);        
    }
    else 
    {
        OxygenStress = 1.;
    }
    
    WatBal->WaterStress = OxygenStress;
    }
