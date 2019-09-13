#include <stdio.h>
#include <math.h>
#include "astro.h"
#include "extern.h"
#include "wofost.h"

#define  ANGLE  -4.0
#define  PI     3.1415926
#define  RAD	0.0174533

// Michaelis-Menten constants for CO2 and O2 at 25oC for C4
//#define  KMC25  650.   // umol/mol
//#define  KMO25  450.   // mmol/mol

// Michaelis-Menten constants for CO2 and O2 at 25oC for C3
# define KMC25  404.9  // umol/mol
# define KMO25  278.4 // mmol/mol
  


// Constants related to the Farquhar-type photosynthesis model
#define O2      210.    // oxygen concentration(mmol/mol)
#define MaxCarboxRate  65330.  // energy of activation for Vcmx(J/mol)
#define EAKMC   79430.  // energy of activation for KMC (J/mol)
#define EAKMO   36380.  // energy of activation for KMO (J/mol)
#define EARD    46390.  // energy of activation for dark respiration(J/mol)
#define DEJMAX  200000. // energy of deactivation for JMAX (J/mol)
#define SJ      650.    // entropy term in JT equation (J/mol/K)
#define PHI2M   0.85    // maximum electron transport efficiency of PS II
#define HH      3.      // number of protons required to synthesise 1 ATP
#define RDVX25  0.0089  // ratio of dark respiration to Vcmax at 25oC

float ScatCoef =0.2;
float XGauss[] ={0.0469101,0.2307534,0.5      ,0.7692465,0.9530899};
float WGauss[] ={0.1184635,0.2393144,0.2844444,0.2393144,0.1184635};


/* ----------------------------------------------------------------------------*/
/*  function InstantAssimilation(float KDiffuse, float EFF, float AssimMax,    */
/*                            float SinB, float PARDiffuse, float PARDirect)   */ 
/*  Purpose: Calculation of the instant Assimilation rate as a function of     */
/*  radiation using the three point Gaussian integration method.               */
/*-----------------------------------------------------------------------------*/

float InstantAssimilation(float KDiffuse, float EFF, float AssimMax, float SinB, 
        float PARDiffuse, float PARDirect)
{
    int i;
    float AbsorbedRadiationDiffuse, AbsorbedRadiationTotal, AbsorbedRadiationDirect;
    float AbsorbedShadedLeaves, AbsorbedDirectLeaves;
    float AssimShadedLeaves, AssimSunlitLeaves, AssimTotal;
    float Reflection, KDirectBl, KDirectTl;
    float GrossCO2, FractionSunlitLeaves, LAIC ;

    /* Extinction coefficients KDIF,KDIRBL,KDIRT */
    Reflection  = (1.-sqrt(1.-ScatCoef))/(1.+sqrt(1.-ScatCoef))*(2/(1+1.6*SinB));
    KDirectBl   = (0.5/SinB)*KDiffuse/(0.8*sqrt(1.-ScatCoef));
    KDirectTl   = KDirectBl*sqrt(1.-ScatCoef);

    /* Three-point Gaussian integration over LAI */
    GrossCO2  = 0.;
    for (i=0;i<3;i++)
    {
       LAIC   = Crop->st.LAI*XGauss[i];
        
       /* Absorbed radiation */
       AbsorbedRadiationDiffuse = (1.-Reflection)*PARDiffuse*KDiffuse * exp(-KDiffuse * LAIC);
       AbsorbedRadiationTotal   = (1.-Reflection)*PARDirect*KDirectTl * exp(-KDirectTl * LAIC);
       AbsorbedRadiationDirect  = (1.-ScatCoef)  *PARDirect*KDirectBl * exp(-KDirectBl * LAIC);

       /* Absorbed flux in W/m2 for shaded leaves and assimilation */
       AbsorbedShadedLeaves = AbsorbedRadiationDiffuse  + AbsorbedRadiationTotal - AbsorbedRadiationDirect;
       AssimShadedLeaves    = AssimMax*(1.-exp (-AbsorbedShadedLeaves*EFF/max(2.0,AssimMax)));

       /* Direct light absorbed by leaves perpendicular on direct */
       /* beam and assimilation of sunlit leaf area               */
       AbsorbedDirectLeaves=(1 - ScatCoef)*PARDirect/SinB;
       if (AbsorbedDirectLeaves <= 0) AssimSunlitLeaves = AssimShadedLeaves;
       else AssimSunlitLeaves = AssimMax*(1. - (AssimMax - AssimShadedLeaves)*
              (1 - exp( -AbsorbedDirectLeaves*EFF/max(2.0,AssimMax)))/(EFF*AbsorbedDirectLeaves));

        /*  Fraction of sunlit leaf area and local assimilation rate  */ 
        FractionSunlitLeaves  = exp(-KDirectBl*LAIC);
        AssimTotal = FractionSunlitLeaves*AssimSunlitLeaves + (1. - FractionSunlitLeaves)*AssimShadedLeaves;

        /*  Integration */
        GrossCO2 += AssimTotal * WGauss[i];
    }
    
    return (GrossCO2 * Crop->st.LAI);     
}

/* ---------------------------------------------------------------------------------*/
/* function ExtinctionDIR()                                                         */
/* Purpose: Calculation of the extinction coefficient for the direct beam radiation */
/* See pg75 formula E2 Crop Systems Dynamics Yin & van Laar 2005                    */
/* ---------------------------------------------------------------------------------*/
float ExtinctionDIR(float *SinB)
{
    float Oav;
    
    // Average projection of the leaves in the solar beam direction
    if (*SinB > sin(Crop->prm.LeafAngle))
    {
        Oav = (*SinB) * cos(Crop->prm.LeafAngle);
    }
    else
    {
        Oav = (*SinB) * cos(Crop->prm.LeafAngle) * 
                 asin(tan(asin(*SinB))/tan(Crop->prm.LeafAngle)) + 
                 sqrt(pow(sin(Crop->prm.LeafAngle,2) - pow((*SinB),2)));
        Oav = Oav * 2./PI;
    }
    
    return Oav/(*SinB);
}


/* ----------------------------------------------------------------------------*/
/* function ExtinctionDIF()                                                    */
/* Purpose: Calculation of the extinction coefficient for diffuse radiation    */
/* See pg75 formula E4 Crop Systems Dynamics Yin & van Laar 2005               */
/* ----------------------------------------------------------------------------*/
float ExtinctionDIF(float *Sc)
{
    float Extinct_15;
    float Extinct_45;
    float Extinct_75;
    
    Extinct_15 = ExtinctionDIR(sin(15.*PI/180.));
    Extinct_45 = ExtinctionDIR(sin(45.*PI/180.));
    Extinct_75 = ExtinctionDIR(sin(75.*PI/180.));
    
    return (-1./Crop->st.LAI) * log(0.178*exp(-Extinct_15 * sqrt(1.-*Sc) * Crop->st.LAI)
                                  + 0.514*exp(-Extinct_45 * sqrt(1.-*Sc) * Crop->st.LAI)
                                  + 0.308*exp(-Extinct_75 * sqrt(1.-*Sc) * Crop->st.LAI));
}


float Reflection(float *K, float *Sc)
{
    float Kbp;
    float Ph;
    /* Scattered beam radiation extinction coefficient */
    Kbp  = (*K) * sqrt(1.-*Sc);

    /* Canopy reflection coefficient for horizontal leaves */
    Ph   = (1. - sqrt(1.- *Sc))/(1. + sqrt(1. - *Sc));

    /* Canopy beam radiation reflection coefficient */
    return 1. - exp(-2. * Ph * (*K)/(1. + (*K)));
}

/* ----------------------------------------------------------------------------*/
/*  function DailyTotalAssimilation()                                          */ 
/*  Purpose: Calculation of the daily assimilation rate using the three point  */
/*  Gaussian integration method.                                               */
/*-----------------------------------------------------------------------------*/
float DailyTotalAssimilation()
{
    int i;
    float KDiffuse, EFF, Factor;
    float Hour, SinB, PAR, NIR, PARDiffuse, PARDirect, AssimMax;
    float NIRDiffuse, NIRDirect;
    float DayTMP;
    float FractionDiffuseRad;
    float K, KdiffPAR, KdiffNIR;
    float DailyTotalAssimilation = 0.;

    KDiffuse = Afgen(Crop->prm.KDiffuseTb, &(Crop->st.Development));

    EFF      = Afgen(Crop->prm.EFFTb, &DayTemp);
    Factor   = Afgen(Crop->prm.CO2EFFTB, &CO2);

    /* Correction for the atmospheric CO2 concentration */
    EFF      = EFF * Factor ;

    AssimMax = Afgen(Crop->prm.FactorAssimRateTemp, &DayTemp) * 
               Afgen(Crop->prm.MaxAssimRate, &(Crop->st.Development)) * 
               Afgen(Crop->prm.CO2AMAXTB, &CO2);

    if (AssimMax > 0. && Crop->st.LAI > 0.)
    {
        for (i=0;i<5;i++)
        {
            Hour       = SunRise + Daylength*XGauss[i];
            SinB       = max (0.,SinLD+CosLD*cos(2.*PI*(Hour+12.)/24.));
            PAR        = 0.5*Radiation[Day][lat][lon]*(SinB*SolarConstant/1367.)/DSinBE;
            NIR        = PAR;
            
            DayTMP = Tmin[Day][lat][lon]+(Tmax[Day][lat][lon]-Tmin[Day][lat][lon]) *
                     sin(PI*(Hour+Daylength/2.-12.)/(Daylength+3.));
            
            /* diffuse light fraction from atmospheric transmission */
            AtmosphTransm = PAR/(0.5*SolarConstant*SinB);
                    
            if (AtmosphTransm > 0.35)
                FractionDiffuseRad = 1.47 - 1.66 * AtmosphTransm;
  
            if (AtmosphTransm <= 0.22 && AtmosphTransm > 0.35)
                FractionDiffuseRad = 1. - 6.4*pow((AtmosphTransm-0.22),2);
  
            if (AtmosphTransm < 0.22)  
                    FractionDiffuseRad = 1.0;
            
            FractionDiffuseRad = max(FractionDiffuseRad, 0.15 + 0.85*(1.-exp(-0.1/SinB)));
            
            /* Incoming diffuse and direct PAR */
            PARDiffuse = PAR * FractionDiffuseRad ;
            PARDirect  = PAR - PARDiffuse;
            
            /* Incoming diffuse and direct NIR */
            NIRDiffuse = NIR * FractionDiffuseRad ;
            NIRDirect  = NIR - NIRDiffuse;
            
            /* Extinction */
            K = ExtinctionDIR(&SinB);
            
            float ScPAR = 0.2;     // leaf scattering coefficient for PAR
            float ScNIR = 0.8;     // leaf scattering coefficient for NIR
            KdiffPAR = ExtinctionDIF(&ScPAR);
            KdiffNIR = ExtinctionDIF(&ScNIR);
            
            CALL REFL (SCPPAR,KB, KBPPAR,PCBPAR)
            CALL REFL (SCPNIR,KB, KBPNIR,PCBNIR)
            
            
            DailyTotalAssimilation = DailyTotalAssimilation + 
                InstantAssimilation(KDiffuse,EFF,AssimMax,SinB,PARDiffuse,PARDirect) * WGauss[i];
        }  
    }
    return(DailyTotalAssimilation*Daylength);
}


/* ----------------------------------------------------------------------------*/
/*  function Correct()                                                         */ 
/*  Purpose: Correct the daily assimilation rate for low temperatures          */
/*-----------------------------------------------------------------------------*/
float Correct(float Assimilation)
{
    int PreviousDay;
    int Counter;
    int number = 7;
    float TminLowAvg = 0.;
    

    if (Crop->GrowthDay < 6)
    {
        number = Crop->GrowthDay;
    }
    
    Counter = 0;
    PreviousDay = Day;
    while (PreviousDay >= 0 && Counter < number)
    {
      TminLowAvg += Tmin[PreviousDay--][lat][lon]; 
      Counter++;
    }

    TminLowAvg = TminLowAvg/Counter;
    return (Assimilation*Afgen(Crop->prm.FactorGrossAssimTemp, &TminLowAvg)*30./44.);

}

float Photosynthese()
{
    float Upar; 
    float TempLeaf;
    float GammaX;
    float Kmc, Kmo;
    float VCT, JT;
    float VCmax, Jmax;
    float FPseud;
    float CO2_leakage;
    float CC, SF, FQ, Fcyc;
    float  Alpha2, X, J2;
    float Vc, Vj;
    float Rdt;
    
    /* Par photonflux umol/m2/s absorbed by leaf photo-sytems */
    Upar = 4.56  * 0.5 *Radiation[Day][lat][lon];
    
    /* Michaelis-Menten constants for CO2 and O2 */
    Kmc = KMC25 *exp((1./298 - 1/(TempLeaf + 273.15)) * EAKMC/8.314);
    Kmo = KMO25 *exp((1./298 - 1/(TempLeaf + 273.15)) * EAKMO/8.314);
    
    //CO2 compensation point in the absence of dark respiration
    GammaX = 0.5*exp(-3.3801+5220./(TempLeaf+273.)/8.314)*O2*Kmc/Kmo;

    //Arrhenius function for the effect of temperature on carboxylation
    VCT    =    exp((1./298.-1./(TempLeaf+273.)) * MaxCarboxRate/8.314);
    
    // function for the effect of temperature on electron transport
    JT = exp((1./298.-1./(TempLeaf + 273.)) * Crop->prm.EnAcJmax/8.314)*
           (1. + exp(SJ/8.314 - Crop->prm.DEJmax/298./8.314))/
           (1. + exp(SJ/8.314-1./(TempLeaf + 273.) *Crop->prm.DEJmax*8.314));
    
    //maximum rates of carboxylation(VCMX) and of electron transport(JMAX)
    VCmax  = Crop->prm.XVN * VCT * ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN);
    Jmax   = Crop->prm.XJN * JT *  ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN);
    
    // CO2 concentration at carboxylation site & electron pathways and
    // their stoichiometries
    FPseud = 0.;           //assuming no pseudocyclic e- transport
    if (!Crop->prm.C3)       //C4
    {
        CO2_leakage   = 0.2;         //CO2 leakage from bundle-sheath to mesophyll
        CC   = 10.*Crop->CO2int;     //to mimic C4 CO2 concentrating mechanism
        SF   = 2.*(CC-GammaX)/(1.-CO2_leakage);
        FQ   = 1.- FPseud- 2.*(4.*CC+8.*GammaX)/HH/(SF+3.*CC+7.*GammaX);
        Fcyc = FQ;
    }
    else
    {
        CC   = Crop->CO2int; 
        SF   = 0.
        FQ   = 0.
        Fcyc = 1.-(FPseud*HH*(SF+3.*CC+7.*GammaX)/(4.*CC + 8.*GammaX)+1.)/
                      (HH*(SF+3.*CC+7.*GammaX)/(4.*CC + 8.*GammaX)-1.);
    }

    // electron transport rate in dependence on PAR photon flux
    Alpha2 = (1.-Fcyc)/(1.+(1.-Fcyc)/PHI2M);
    X      = Alpha2 * Upar/max(1.E-10,Jmax);
    J2     = Jmax*(1+X-((1+X)**2-4.*X*Crop->prm.Theta)**0.5)/2./Crop->prm.Theta;

    // rates of carboxylation limited by Rubisco and electron transport
    Vc   = VCmax * CC/(CC + Kmc*(O2/Kmo+1.));
    Vj   = J2 * CC*(2.+FQ-Fcyc)/HH/(SF+3.*CC+7.*GammaX)/(1.-Fcyc);

    // gross rate of leaf photosynthesis
    Crop->rt.LeafPhoto  = max(1.E-10, (1.E-6)*44. * (1.-GammaX/CC)*min(Vc,Vj));

    // leaf dark respiration rate
    RDVX25 = 0.0089;      //ratio of dark respiration to Vcmax at 25oC
    Rdt    = exp((1./298.-1./(TempLeaf + 273.))*EARD/8.314);
    Crop->rt.DarkResp = (1.E-6)*44. * RDVX25 * (Crop->prm.XVN * 
            ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN)) * Rdt;    
    
}