#include <stdio.h>
#include <stdlib.h>
#include "wofost.h"
#include "extern.h"

/* ----------------------------------------------------------------------------------*/
/*  function Conversion                                                              */
/*  Purpose: Convert the net assimilation products into plant dry matter kg ha-1 d-1 */
/* ----------------------------------------------------------------------------------*/

float Conversion(float NetAssimilation)
{
    float fr, root, shoots;

    fr    = Afgen(Crop->prm.Roots, &(Crop->st.Development));
    root  = fr/Crop->prm.ConversionRoots;
    shoots =  Afgen(Crop->prm.Stems, &(Crop->st.Development))/Crop->prm.ConversionStems;
    shoots += Afgen(Crop->prm.Leaves, &(Crop->st.Development))/Crop->prm.ConversionLeaves;	
    shoots += Afgen(Crop->prm.Storage, &(Crop->st.Development))/Crop->prm.ConversionStorage;

    /* conversion */
    return NetAssimilation/(shoots*(1-fr)+root);
}