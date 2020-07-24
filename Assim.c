#include <stdio.h>
#include <math.h>
#include "astro.h"
#include "extern.h"
#include "wofost.h"
#include "penman.h"

#define  ANGLE  -4.0
#define  PI     3.1415926
#define  RAD	0.0174533

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

#define ScatPAR   0.2       // leaf scattering coefficient for PAR
#define ScatNIR   0.8       // leaf scattering coefficient for NIR  
#define PcdPAR  0.057    // canopy diffuse PAR reflection coefficient
#define PcdNIR  0.389    // canopy diffuse NIR reflection coefficient  

//float ScatCoef =0.2;
float XGauss[] ={0.0469101,0.2307534,0.5      ,0.7692465,0.9530899};
float WGauss[] ={0.1184635,0.2393144,0.2844444,0.2393144,0.1184635};

void PhotoAciveN(float Kn, float KDiffuse, float Slnt,float *NP_su, float *NP_sh)
{
    float NP_c;
    // Total photosynthetic nitrogen in canopy
    NP_c   = Slnt*(1.-exp(-Kn * Crop->st.LAI))/Kn- Crop->prm.SLMIN*Crop->st.LAI;

    // Photosynthetic nitrogen for sunlit and shaded parts of canopy
    *NP_su  = Slnt*(1.-exp(-(Kn+KDiffuse)*Crop->st.LAI))/(Kn+KDiffuse) - 
            Crop->prm.SLMIN*(1.-exp(-KDiffuse*Crop->st.LAI))/KDiffuse;
    *NP_sh  = NP_c - NP_su;   
}

float InternalCO2(float *Temperature)
{
    float SatVapLeaf, VapDeficitLv;
    float Kmc;
    float Kmo;
    float Gamma, Gamma0, GammaX;
    float RatioDarkResVCX;
    float FVPD; //Slope for linear effect of VPD on CO2internal/CO2ambient
    float RatioCO2intCO2amb;
    
    FVPD   = insw(Crop->C3C4, 0.195127, 0.116214);
    
    // Air-to-leaf vapour pressure deficit
    SatVapLeaf   = 0.611 * exp(17.4 * *Temperature / (*Temperature + 239.));
    VapDeficitLv = max(0., SatVapLeaf - Vapour[Day][lat][lon]);

    Kmc    = KMC25*exp((1./298.-1./(*Temperature + 273.))*EAKMC/8.314);
    Kmo    = KMO25*exp((1./298.-1./(*Temperature + 273.))*EAKMO/8.314);
    GammaX = 0.5*exp(-3.3801+5220./(*Temperature + 273.)/8.314)*O2*Kmc/Kmo;

    // CO2 compensation point (GAMMA)
    RatioDarkResVCX  = RDVX25 * exp((1./298.-1./(Crop->st.LeafTemp+273.))*(EARD-MaxCarboxRate)/8.314);
    Gamma0 = (GammaX + RatioDarkResVCX * Kmc * (1.+O2/Kmo))/(1. - RatioDarkResVCX);
    Gamma  = insw(Crop->C3C4, Gamma0/10., Gamma0);

    // Internal/ambient CO2 ratio, based on data of Morison & Gifford (1983)
    RatioCO2intCO2amb  = 1.-(1.-Gamma/CO2)*(0.14 + FVPD * VapDeficitLv);

    // Intercellular CO2 concentration
    return RatioCO2intCO2amb * CO2;
}

float LeafPhotoResp(float *PAR, float *Temperature, float *LeafPhoto, float *DarkResp)
{
    float Upar; 
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
    Upar = 4.56  * *PAR;
    
    /* Michaelis-Menten constants for CO2 and O2 */
    Kmc = KMC25 *exp((1./298 - 1/(*Temperature + 273.15)) * EAKMC/8.314);
    Kmo = KMO25 *exp((1./298 - 1/(*Temperature + 273.15)) * EAKMO/8.314);
    
    //CO2 compensation point in the absence of dark respiration
    GammaX = 0.5*exp(-3.3801+5220./(*Temperature+273.)/8.314)*O2*Kmc/Kmo;

    //Arrhenius function for the effect of temperature on carboxylation
    VCT    =    exp((1./298.-1./(*Temperature+273.)) * MaxCarboxRate/8.314);
    
    // function for the effect of temperature on electron transport
    JT = exp((1./298.-1./(*Temperature + 273.)) * Crop->prm.EnAcJmax/8.314)*
           (1. + exp(SJ/8.314 - Crop->prm.DEJmax/298./8.314))/
           (1. + exp(SJ/8.314-1./(*Temperature + 273.) *Crop->prm.DEJmax*8.314));
    
    //maximum rates of carboxylation(VCMX) and of electron transport(JMAX)
    VCmax  = Crop->prm.XVN * VCT * ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN);
    Jmax   = Crop->prm.XJN * JT *  ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN);
    
    // CO2 concentration at carboxylation site & electron pathways and
    // their stoichiometries
    FPseud = 0.;             //assuming no pseudocyclic e- transport
    if (Crop->C3C4 == -1)    //C4
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
    *LeafPhoto  = max(1.E-10, (1.E-6)*44. * (1.-GammaX/CC)*min(Vc,Vj));

    // leaf dark respiration rate
    RDVX25 = 0.0089;      //ratio of dark respiration to Vcmax at 25oC
    Rdt    = exp((1./298.-1./(*Temperature + 273.))*EARD/8.314);
    *DarkResp = (1.E-6)*44. * RDVX25 * (Crop->prm.XVN * 
            ((Crop->N_st.leaves/Crop->st.LAI) - Crop->prm.SLMIN)) * Rdt;    
}


float Photosynthesis(float Rad, float Temp)
{
    float ICO2;
    float AirLeafVPD;
    float CndCO2;
    float VPD, Vap, Svp, Svp_L, Delta, Delta_L, Dif, TempLeaf;
    float LeafPhoto, DarkResp, RnetAbs;
    float PT;
    
    float Lhvap  = 2.45e6;  // Latent heat of evaporation of water (J/kg=J/mm)  
    float Vhca   = 1200.;   // Volumetric heat capacity
    
    LeafPhoto = 0.;
    DarkResp  = 0.;
    RnetAbs   = 0.;
    ICO2 = InternalCO2(&Temp);
    
    Vap       = 0.1 * (Vapour[Day][lat][lon]);
    // Air to leaf vapour pressure deficit
    Svp       =  0.6108 * exp((17.27 * Temp) / (237.3 + Temp)); //eq 8 pg 13
    VPD       = max(0.,Svp - Vap);
    Delta     = (4098. * Svp)/pow((Temp + 237.3), 2);
    
    LeafPhotoResp(&Rad, &Temp, &LeafPhoto, &DarkResp);
    
    // Potential conductance for CO2
    CndCO2 = (LeafPhoto - DarkResp)*((273.+Temp)/0.53717)/(CO2 - ICO2); //eq 4 pg11
    
    // potential stomatal resistance to water
    R_stomata = max(1e-10, 1./CndCO2 - RB_Water_sh *1.3 - RT_canopy)/1.6;
    
    CalcPenmanMonteith(&Rad, &RnetAbs, float *PT, float *Frac, &R_stomata);
    
    Dif = limit(-25.,25., RnetAbs - Lhvap * PT) * (RT_canopy + RB_Heat);
    TempLeaf = Temp +Dif;
    
    // Second round to determine the final photosynthesis and transpiration
    Svp_L = 0.6108 * exp((17.27 * TempLeaf) / (237.3 + Temp));
    ICO2  = InternalCO2(TempLeaf);
    
    // Potential conductance for CO2
    CndCO2 = (LeafPhoto - DarkResp)*((273.+TempLeaf)/0.53717)/(CO2 - ICO2); //eq 4 pg11
    
    // potential stomatal resistance to water
    R_stomata = max(1e-10, 1./CndCO2 - RB_Water_sh *1.3 - RT_canopy)/1.6;
    LeafPhotoResp(Rad, TempLeaf);
    
    if (TempLeaf == Temp)
        Delta_L = 1.;
    else
        Delta_L = (Svp_L - Svp)/(TempLeaf - Temp);    //eq 6 pg 12
    
    CalcPenmanMonteith(&R_stomata);
    
    
    
}

/* ----------------------------------------------------------------------------*/
/*  function InstantAssimilation(float KDiffuse, float EFF, float AssimMax,    */
/*                            float SinB, float PARDiffuse, float PARDirect)   */ 
/*  Purpose: Calculation of the instant Assimilation rate as a function of     */
/*  radiation using the three point Gaussian integration method.               */
/*-----------------------------------------------------------------------------*/

float InstantAssimilation()
{
    int i;

    

     
    (AssimShaded + AssimSunlit);
    
    return (AssimShaded + AssimSunlit);     
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
    float Hour, SinB, PAR, PARDiffuse, PARDirect, AssimMax; 
    float NIR, NIRDirect, NIRDiffuse;
    float DailyTotalAssimilation = 0.;
    float Tx, Tn;
    float GrossCO2;
    
    float ReflectPAR, KDirBlPAR, KDirTlPAR;
    float ReflectNIR, KDirBlNIR, KDirTlNIR;

    float Kwind, Kn; // Extinction coefficient wind and nitrogen
    float Kln, Nbk;
    float Sln, Slnt, Slnnt;
    float PARsoil, NIRsoil;
    
    float BndConducH;   
    float BndConducCnp;  
    float BndConducSunlit;
    float BndConducShaded;
    float KDifPAR, KDifNIR; 
    
    GrossCO2  = 0.;
    
    KDiffuse = Afgen(Crop->prm.KDiffuseTb, &(Crop->st.Development));
    
    if (Crop->st.LAI > 0.)
    {
        for (i=0;i<5;i++)
        {
            Hour       = 12.0 + 0.5 * Daylength * XGauss[i];
            SinB       = max (0.,SinLD+CosLD*cos(2.*PI*(Hour+12.)/24.));
            Tn         = Tmin[Day][lat][lon];
            Tx         = Tmax[Day][lat][lon];
            Temp       = Tn + (Tx-Tn) * sin(PI*(Hour + Daylength/2.-12.)/(Daylength + 3.));          
            PAR        = 0.5*Radiation[Day][lat][lon]*SinB*(1.+0.4*SinB)/DSinBE;
            PARDiffuse = min (PAR,SinB*DiffRadPP);
            PARDirect  = PAR - PARDiffuse;
            NIR        = PAR;
            NIRDiffuse = PARDiffuse;
            NIRDirect  = NIR - NIRDiffuse;
            
            // Taken from the original wofost code and 
            // Spitters 1986: Agricultural and Forest Meteorology, 38 (1986) 231--242
            // Extinction coefficients KDIF,KDIRBL,KDIRT PAR 
            KDifPAR     = 0.8*sqrt(1-ScatPAR);
            ReflectPAR  = (1.-sqrt(1.-ScatPAR))/(1.+sqrt(1.-ScatPAR))*(2/(1+1.6*SinB));
            KDirBlPAR   = (0.5/SinB)*KDiffuse/KDifPAR;
            KDirTlPAR   = KDirBlPAR*sqrt(1.-ScatPAR);    // PAR extinction coefficient of black leaves assuming spherical angle distribution


            /* Extinction coefficients KDIF,KDIRBL,KDIRT NIR */
            ReflectNIR  = (1.-sqrt(1.-ScatNIR))/(1.+sqrt(1.-ScatNIR))*(2/(1+1.6*SinB));
            KDifNIR     = 0.8*sqrt(1-ScatNIR);
            KDirBlNIR   = (0.5/SinB)*KDiffuse/KDifNIR;
            KDirTlNIR   = KDirBlPAR*sqrt(1.-ScatNIR);    // NIR extinction coefficient of black leaves

            // Assume extinction coefficient wind = extinction coefficient PAR
            Kwind = KDiffuse;

            // Extinction coefficient of nitrogen  */
            Kln    = KDiffuse * (Crop->N_st.leaves - Crop->prm.SLMIN * Crop->st.TLAI);
            Nbk    = Crop->prm.SLMIN * (1.-exp(-KDiffuse * Crop->st.LAI));
            Kn     = 1./Crop->st.TLAI*log((Kln + Nbk)/(Kln*exp(-KDiffuse * Crop->st.TLAI)+Nbk)); // eq41 pg35

            // Specific leaf nitrogen and its profile in the canopy
            Sln   = Crop->N_st.leaves/Crop->st.LAI;
            Slnt  = Crop->N_st.leaves * Kn;
            Slnnt = (Crop->N_st.leaves + 0.001*Crop->N_st.leaves)*Kn /(1.-exp(-Kn*Crop->st.LAI));   

            // Turbulence resistance for canopy (RT) and for soil (RTS) 
            RT_canopy  = 0.74*pow((log((2.-0.7*Crop->st.Height)/(0.1*Crop->st.Height))),2)/ (0.16 * Windspeed[Day][lat][lon]);
            RT_soil    = 0.74*pow((log(56.)),2)/(0.16 * Windspeed[Day][lat][lon])

            /* Canopy boundary layer conductance for sunlit and shaded leaves */
            BndConducH       = 0.01*sqrt(Windspeed[Day][lat][lon]/Crop->prm.LeafWidth);
            BndConducCnp     = (1.-exp(- 0.5 * Kwind * Crop->st.LAI))/(0.5 * Kwind ) * BndConducH;
            BndConducSunlit  = (1.-exp(-(0.5 * Kwind + KDirTlPAR) * Crop->st.LAI)) /
                    (0.5 * Kwind + KDirTlPAR) * BndConducH; 
            BndConducShaded  = BndConducCnp - BndConducSunlit;

            // Canopy boundary layer resistance for sunlit and shaded leaves 
            RB_Heat_su  = 1./BndConducSunlit;       // boundary layer resistance to heat,sunlit part
            RB_Water_su = 0.93 * RB_Heat_su;        // boundary layer resistance to H2O, sunlit part
            RB_Heat_sh  = 1./BndConducShaded;       // boundary layer resistance to heat,shaded part
            RB_Water_sh = 0.93 * RB_Heat_sh;        // boundary layer resistance to H2O, shaded part 

            // Boundary layer resistance for soil
            RS_Heat   = 172.*sqrt(0.05/max(0.1,Windspeed[Day][lat][lon]*exp(-Kwind * Crop->st.TLAI)));
            RS_Water  = 0.93*RS_Heat

            //  Fraction of sunlit leaf area 
            FracSunlit = (1./(KDirBlPAR * Crop->st.LAI)) * (1 - exp(-KDirBlPAR*Crop->st.LAI));
            FracShaded = 1 - FracSunlit;

            // Absorbed PAR radiation */
            APARDiffuse = (1.-ReflectPAR)*PARDiffuse*KDiffuse * exp(-KDiffuse * Crop->st.LAI);
            APARTotal   = (1.-ReflectPAR)*PARDirect*KDirTlPAR * exp(-KDirTlPAR * Crop->st.LAI);
            APARDirect  = (1.-ScatPAR)   *PARDirect*KDirBlPAR * exp(-KDirBlPAR * Crop->st.LAI);

            // Absorbed NIR radiation */
            ANIRDiffuse = (1.-ReflectNIR)*NIRDiffuse*KDiffuse * exp(-KDiffuse * Crop->st.LAI);
            ANIRTotal   = (1.-ReflectNIR)*NIRDirect*KDirTlNIR * exp(-KDirTlNIR * Crop->st.LAI);
            ANIRDirect  = (1.-ScatNIR)   *NIRDirect*KDirBlNIR * exp(-KDirBlNIR * Crop->st.LAI);

            // Absorbed flux in W/m2 for shaded leaves  */
            APARShadedLeaves = APARDiffuse  + APARTotal - APARDirect;
            ANIRShadedLeaves = ANIRDiffuse  + ANIRTotal - ANIRDirect;
            AShadedLeaves    = APARShadedLeaves + ANIRShadedLeaves

            // Direct light absorbed by leaves perpendicular on direct 
            // beam and assimilation of sunlit leaf area               
            APARSunlitLeaves = (1 - ScatPAR) * PARDirect/SinB;
            ANIRSunlitLeaves = (1 - ScatNIR) * NIRDirect/SinB;
            ASunlitLeaves = APARSunlitLeaves + ANIRSunlitLeaves;

             // Fraction of sunlit leaf area and local assimilation rate  
             // FractionSunlitLeaves  = exp(-KDirectBl*LAIC);
             AssimShaded = Photosynthesis(AShadedLeaves, Temp) 
             AssimSunlit = Photosynthesis(ASunlitLeaves, Temp);

            // Absorbed total radiation (PAR+NIR) by soil
            PARsoil  = 0.1;                                   //soil PAR reflection
            if (WatBal->st.Moisture - 0.5 < 0.)
                NIRsoil =  0.52-0.68 * WatBal->st.Moisture;   //soil NIR reflection
            else
                NIRsoil = 0.18;

            //  Absorbed radiation soil
            ASoil = (1.-PARsoil)*(PARDirect * exp(-ReflectPAR * Crop->st.TLAI) + 
                     PARDiffuse * exp(-KDifPAR * Crop->st.TLAI)) +
                     (1.-NIRsoil)*(NIRDirect * exp(-ReflectNIR * Crop->st.TLAI) +
                     NIRDiffuse *exp(-KDifNIR * Crop->st.TLAI));

            InstantAssimilation(FracSunlit,APARSunlitLeaves);
            InstantAssimilation(FracShaded,APARShadedLeaves);
            
            
            DailyTotalAssimilation = DailyTotalAssimilation + 
                InstantAssimilation(KDiffuse, Temp, SinB, PARDiffuse, PARDirect, NIRDiffuse, NIRDirect) * WGauss[i];
        }  
    }
    return(DailyTotalAssimilation*Daylength);
}


