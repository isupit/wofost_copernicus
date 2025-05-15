#include <stdio.h>
#include <stdlib.h>
#include "extern.h"
#include "wofost.h"

void header(FILE *fp)
{
    //fprintf(fp,"Lat   Lon    Sowing length TSM1 TSM2  Avg    Adev   Stdev  Var      Skew  Curt\n");
    fprintf(fp,"Lat,Lon,Sowing,Length,TSM1,TSM2,Avg,Adev,Stdev,Var,Skew,Curt\n");
}

void Output(FILE *fp ) 
{
    float ave, adev, sdev, var, skew, curt, lngth;
    int mnth, dy, dekad, i;
    
    sscanf(Grid->start,"%2d-%2d",&mnth,&dy);
   
    if (mnth < 1 || mnth > 12) exit(0);
    if (dy < 1 || dy >31) exit (0);
     
    // convert sowing day to dekad
    dekad = (mnth -1) * 3;
    if (dy <= 10)
        dekad += 1;
    else if(dy <=20)
        dekad += 2;
    else
        dekad += 3;   
                  
    
    if (Crop->Seasons > 2) {
        lngth =0;
        for (i= 1; i <= Crop->Seasons; i++) {
            lngth = lngth + Grid->length[i];
    }
       
        lngth = lngth/Crop->Seasons;
               
        Moment(Grid->twso, Crop->Seasons, &ave, &adev, &sdev, &var, &skew, &curt);  
       
        fprintf(fp, "%4.2f, %4.2f, %6d, %6.0f, %4.0f, %4.0f, %6.0f, %6.0f, %6.0f, %9.0f, %5.2f, %5.2f %5d\n", 
        Latitude[Lat], 
        Longitude[Lon], 
        dekad,   //dekad sowing
        lngth,
        Crop->prm.TempSum1,
        Crop->prm.TempSum2,
        ave,
        adev,
        sdev,
        var,
        skew,
        curt,
        Crop->Seasons);
    }
 }

    
//void Output(FILE *fp)
//{    
//    fprintf(fp,"%7.2f\t%7.2f\t\t%4d\t\t%3d\t\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t"
//            "\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n",
//        Latitude[Lat],
//        Longitude[Lon],
//        MeteoYear[Day],
//        MeteoDay[Day],
//        Crop->st.Development,
//        Crop->st.leaves,
//        Crop->st.stems,
//        Crop->st.storage,
//        Crop->st.roots,
//        Crop->st.LAI,
//        WatBal->WaterStress,
//        WatBal->st.Moisture,
//        Rain[Lon][Lat][Day],
//        WatBal->rt.Infiltration,
//        WatBal->rt.Runoff,
//        WatBal->rt.Loss,
//        Crop->N_st.Indx,
//        Crop->P_st.Indx,
//        Crop->K_st.Indx,
//        Site->st_N_tot,
//        Site->st_P_tot,
//        Site->st_K_tot,
//        Crop->N_rt.Uptake,
//        Crop->P_rt.Uptake,
//        Crop->K_rt.Uptake,
//        Crop->N_rt.Demand_lv + Crop->N_rt.Demand_st + Crop->N_rt.Demand_ro,
//        Crop->P_rt.Demand_lv + Crop->P_rt.Demand_st + Crop->P_rt.Demand_ro,
//        Crop->K_rt.Demand_lv + Crop->K_rt.Demand_st + Crop->K_rt.Demand_ro);
//}

