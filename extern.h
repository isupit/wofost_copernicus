/* General global parameters */
#include "wofost.h"

#ifndef EXTERN_H
#define EXTERN_H

float Temp;
float DayTemp;

/* General help functions */
extern float max(float a, float b);
extern float min(float a, float b);
extern float Afgen();
extern float List();
extern float limit(float a, float b, float c);
extern float notnul(float x);
extern float insw(float x1, float x2, float x3);
extern int leap_year(int year);
extern void RatesToZero();
extern void Output();
extern void header();
extern int GetSimInput();
extern void GetMeteoInput();
extern void IfSowing();

/* Additional functions */
extern int Astro();
extern void CalcPenman();
extern void CalcPenmanMonteith();
extern void Clean();

/* Crop growth */
extern void Partioning();
extern void HeatStress();
extern void RateCalculationCrop();
extern void Growth(float NewPlantMaterial);
extern void IntegrationCrop();
extern void InitializeCrop();
extern int EmergenceCrop(int Emergence);

extern void DevelopmentRate();
extern void LeaveGrowth();
extern float DailyTotalAssimilation();
extern float DyingLeaves();
extern float InstantAssimilation(float KDiffuse, float EFF, float AssimMax, float SinB, float PARDiffuse, float PARDirect);
extern float LeaveAreaIndex();
extern float Correct(float GrossAssimilation);
extern float RespirationRef(float TotalAssimilation);
extern float Conversion(float NetAssimilation);


/* Nutrients */
extern void CropNutrientRates();
extern void InitializeNutrients();
extern void IntegrationNutrients();
extern void NutritionINDX();
extern void NutrientLoss();
extern void NutrientMax();
extern void NutrientPartioning();
extern void NutrientRates();
extern void NutrientOptimum();
extern void NutrientDemand();
extern void SoilNutrientRates();
extern void NutrientTranslocation();
extern void RateCalcultionNutrients();



/* Read data */
extern int GetAfgenTables();
extern int GetMeteoData();
extern int FillAfgenTables();

extern void GetCropData();
extern void FillCropVariables();
extern void GetSiteData();
extern void FillSiteVariables();
extern void  GetSoilData();
extern void FillSoilVariables();
extern void GetManagement();
extern void FillManageVariables();

/* Water balance */
extern void InitializeWatBal();
extern void RateCalulationWatBal();
extern void IntegrationWatBal();
extern void EvapTra();

#endif	// EXTERN_H

