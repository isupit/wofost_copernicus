#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wofost.h"
#include "extern.h"

/* -----------------------------------------------------------------------------------------*/
/*  function Partioning()                                                                   */
/*  Purpose: Calculate the partioning factors and correct them for nutrient or water stress */ 
/* -----------------------------------------------------------------------------------------*/

void Partioning()
{
    float factor;
       
    factor = fmax(1., 1./(WatBal->WaterStress + 0.5));
    Crop->fac_ro = fmin(0.6, Afgen(Crop->prm.Roots, &(Crop->st.Development)) * factor);
    Crop->fac_lv = Afgen(Crop->prm.Leaves, &(Crop->st.Development));
    Crop->fac_st = Afgen(Crop->prm.Stems, &(Crop->st.Development));
    Crop->fac_so = Afgen(Crop->prm.Storage, &(Crop->st.Development));
}	
