#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astro.h"
#include "extern.h"
#include "wofost.h"
#include "penman.h"
#include "assim.h"

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

#define ScatPAR    0.2      // leaf scattering coefficient for PAR
#define ScatNIR    0.8      // leaf scattering coefficient for NIR  
#define ReflPAR_d  0.057    // canopy diffuse PAR reflection coefficient
#define ReflNIR_d  0.389    // canopy diffuse NIR reflection coefficient 


//float ScatCoef =0.2;
float XGauss[] ={0.0469101,0.2307534,0.5      ,0.7692465,0.9530899};
float WGauss[] ={0.1184635,0.2393144,0.2844444,0.2393144,0.1184635};

float KBeam(float sin_b)
{
    float B, OAV, LeafAngle;
    
    LeafAngle = Crop->prm.LeafAngle * RAD;
    
    // Solar elevation in radians
     B      = asin(sin_b);

    // Average projection of leaves in the direction of a solar beam
    if (sin_b >= sin(LeafAngle))
        OAV = sin_b * cos(LeafAngle);
    else
        OAV = 2./3.141592654*(sin_b * cos(LeafAngle) * asin(tan(B)/tan(LeafAngle))
             + sqrt(pow(sin(LeafAngle),2) - pow(sin_b,2)));
    
    // Beam radiation extinction coefficient
    return(OAV/SinB);
}
   
float KDiff(float Scat)
{    
    float Kb15, Kb45, Kb75;
    
    // Extinction coefficient of beam lights from 15, 45 and 75o elevations
    Kb15 = KBeam(15. * RAD);
    Kb45 = KBeam(45. * RAD);
    Kb75 = KBeam(75. * RAD);

    // Diffuse light extinction coefficient
    return(-1./Crop->st.TLAI*log(0.178*exp(-Kb15*sqrt(1.-Scat)*Crop->st.TLAI)
                    + 0.514*exp(-Kb45*sqrt(1.-Scat)*Crop->st.TLAI)
                    + 0.308*exp(-Kb75*sqrt(1.-Scat)*Crop->st.TLAI)));
}

float Reflect(float SCP)
{ 
    float Ph;
    // Canopy reflection coefficient for horizontal leaves
    Ph  = (1.-sqrt(1.-SCP))/(1.+sqrt(1.-SCP));

    // Canopy beam radiation reflection coefficient
    return (1.-exp(-2.*Ph*Kb/(1.+Kb)));
}

void PhotoAciveN()
{
    float NP_tot, NP_tot_n;
    
    // Total photosynthetic nitrogen in canopy
    NP_tot   = ( Slnt*(1.-exp(-KNitro * Crop->st.LAI))/KNitro - 
            Crop->prm.SLMIN*Crop->st.LAI); 
    NP_tot_n = ( Slnnt*(1.-exp(-KNitro * Crop->st.LAI))/KNitro - 
            Crop->prm.SLMIN*Crop->st.LAI);
    
   // Photosynthetic nitrogen for sunlit and shaded parts of canopy
    Su->NP   = (Slnt*(1.-exp(-(KNitro + Kb) * Crop->st.LAI))/(KNitro + Kb) - 
            Crop->prm.SLMIN*(1.-exp(-Kb * Crop->st.LAI))/Kb);
    Su->NP_n = (Slnnt*(1.-exp(-(KNitro + Kb) * Crop->st.LAI))/(KNitro + Kb) - 
            Crop->prm.SLMIN*(1.-exp(-Kb * Crop->st.LAI))/Kb);
    
    Sh->NP   = NP_tot - Su->NP;
    Sh->NP_n = NP_tot_n - Su->NP_n;
}

void AbsorbedLight_tot()
{
    float APAR_tot, ANIR_tot;
    
    APAR_tot = (1. - ReflPAR_b) * PARDirect *(1.-exp(-KbPAR*Crop->st.LAI)) +
               (1. - ReflPAR_d) * PARDiffuse*(1.-exp(-KdpPAR*Crop->st.LAI));
    
    ANIR_tot = (1. - ReflNIR_b) * NIRDirect *(1.-exp(-KbNIR*Crop->st.LAI)) +
               (1. - ReflNIR_d) * NIRDiffuse*(1.-exp(-KdpNIR*Crop->st.LAI));
    
    Su->APAR = (1.-ScatPAR)  * PARDirect * (1.-exp(-Kb * Crop->st.LAI)) +
              (1.-ReflPAR_d) * PARDiffuse/(KdpPAR+Kb) * KdpPAR*(1.-exp(-(KdpPAR+Kb)*Crop->st.LAI)) + 
              PARDirect * ((1.-ReflPAR_b)/(KbPAR+Kb) * KbPAR *
              (1.-exp(-(KbPAR+Kb) * Crop->st.LAI))-(1.-ScatPAR)*(1.-exp(-2.*Kb*Crop->st.LAI))/2.);
    
    Su->ANIR = (1.-ScatNIR) * NIRDirect * (1.-exp(-Kb * Crop->st.LAI)) +
              (1.-ReflNIR_d) * NIRDiffuse/(KdpNIR+Kb) * KdpNIR*(1.-exp(-(KdpNIR+Kb)*Crop->st.LAI)) + 
              NIRDirect * ((1.-ReflNIR_b)/(KbNIR+Kb) * KbNIR *
              (1.-exp(-(KbNIR+Kb) * Crop->st.LAI))-(1.-ScatNIR)*(1.-exp(-2.*Kb*Crop->st.LAI))/2.);
    
    Su->APAR = APAR_tot - Su->APAR;
    Sh->ANIR = ANIR_tot - Su->ANIR;        
}

void InternalCO2()
{
    float Kmc;
    float Kmo;
    float Gamma, Gamma0, GammaX;
    float RatioDarkResVCX;
    float RatioCO2intCO2amb;
    float Temperature;
    
    Temperature = DTemp + Dif;
    
    // Air-to-leaf vapour pressure deficit
    SatVap     = 0.611 * exp(17.4 * Temperature / (Temperature + 239.));
    VapDeficit = max(0., SatVap - Vapour[Day][lat][lon]);
    Delta      = (4098. * SatVap)/pow((Temperature + 237.3), 2);

    Kmc    = KMC25*exp((1./298.-1./(Temperature + 273.))*EAKMC/8.314);
    Kmo    = KMO25*exp((1./298.-1./(Temperature + 273.))*EAKMO/8.314);
    GammaX = 0.5*exp(-3.3801+5220./(Temperature + 273.)/8.314)*O2*Kmc/Kmo;

    // CO2 compensation point (GAMMA)
    RatioDarkResVCX  = RDVX25 * exp((1./298.-1./(Temperature + 273.))*(EARD-MaxCarboxRate)/8.314);
    Gamma0 = (GammaX + RatioDarkResVCX * Kmc * (1.+O2/Kmo))/(1. - RatioDarkResVCX);
    Gamma  = insw(Crop->C3C4, Gamma0/10., Gamma0);

    // Internal/ambient CO2 ratio, based on data of Morison & Gifford (1983)
    RatioCO2intCO2amb  = 1.-(1.-Gamma/CO2)*(0.14 + Fvpd * VapDeficit);

    // Intercellular CO2 concentration
    CO2int = RatioCO2intCO2amb * CO2;
}

void LeafPhotoResp(float APAR, float NP, float *LeafPhoto, float *DarkResp)
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
    float Temperature;
    
    Temperature = DTemp + Dif;
    
    // Par photonflux umol/m2/s absorbed by leaf photo-sytems 
    Upar = 4.56  * APAR;
    
    // Michaelis-Menten constants for CO2 and O2 */
    Kmc = KMC25 *exp((1./298 - 1/(Temperature + 273.15)) * EAKMC/8.314);
    Kmo = KMO25 *exp((1./298 - 1/(Temperature + 273.15)) * EAKMO/8.314);
    
    // CO2 compensation point in the absence of dark respiration
    GammaX = 0.5*exp(-3.3801+5220./(Temperature +273.)/8.314)*O2*Kmc/Kmo;

    // Arrhenius function for the effect of temperature on carboxylation
    VCT    =    exp((1./298.-1./(Temperature + 273.)) * MaxCarboxRate/8.314);
    
    // Function for the effect of temperature on electron transport
    JT = exp((1./298.-1./(Temperature + 273.)) * Crop->prm.EnAcJmax/8.314)*
           (1. + exp(SJ/8.314 - Crop->prm.DEJmax/298./8.314))/
           (1. + exp(SJ/8.314-1./(Temperature + 273.) *Crop->prm.DEJmax*8.314));
    
    // Maximum rates of carboxylation(VCMX) and of electron transport(JMAX)
    VCmax  = Crop->prm.XVN * VCT * NP;
    Jmax   = Crop->prm.XJN * JT *  NP;
    
    // CO2 concentration at carboxylation site & electron pathways and
    // their stoichiometries
    FPseud = 0.;             //assuming no pseudocyclic e- transport
    if (Crop->C3C4 == -1)    //C4
    {
        CO2_leakage   = 0.2;     //CO2 leakage from bundle-sheath to mesophyll
        CC   = 10. * CO2int;     //to mimic C4 CO2 concentrating mechanism
        SF   = 2.*(CC - GammaX)/(1. - CO2_leakage);
        FQ   = 1.- FPseud- 2.*(4.*CC+8.*GammaX)/HH/(SF+3.*CC+7.*GammaX);
        Fcyc = FQ;
    }
    else
    {
        CC   = CO2int; 
        SF   = 0.;
        FQ   = 0.;
        Fcyc = 1.-(FPseud*HH*(SF+3.*CC+7.*GammaX)/(4.*CC + 8.*GammaX)+1.)/
                      (HH*(SF+3.*CC+7.*GammaX)/(4.*CC + 8.*GammaX)-1.);
    }

    // Electron transport rate in dependence on PAR photon flux
    Alpha2 = (1.-Fcyc)/(1.+(1.-Fcyc)/PHI2M);
    X      = Alpha2 * Upar/max(1.E-10,Jmax);
    J2     = Jmax*(1+X -sqrt( pow((1+X),2) - 4.*X*Crop->prm.Theta))/2./Crop->prm.Theta;

    // Rates of carboxylation limited by Rubisco and electron transport
    Vc   = VCmax * CC/(CC + Kmc*(O2/Kmo+1.));
    Vj   = J2 * CC*(2.+FQ-Fcyc)/HH/(SF+3.*CC+7.*GammaX)/(1.-Fcyc);

    // Gross rate of leaf photosynthesis
    *LeafPhoto  = max(1.e-10, (1.e-6)*44. * (1.-GammaX/CC)*min(Vc,Vj));

    // leaf dark respiration rate
    Rdt    = exp((1./298.-1./(Temperature + 273.))*EARD/8.314);
    *DarkResp = (1.e-6)*44. * RDVX25 * (Crop->prm.XVN * NP) * Rdt;    
}


void InstantAssimTransp(SUSH S)
{
    float R_turb;
    float CndCO2;
    float Svp, Svp_L;
    
    S->R_stomata = 0.;
    R_turb = S->Frac * S->RT;
          
    InternalCO2();
    Svp  = SatVap;
    
    S->ARAD = S->APAR + S->ANIR; 
    LeafPhotoResp(S->ARAD, S->NP, &S->LeafPhoto, &S->DarkResp_n);
    
    // Potential conductance for CO2
    CndCO2 = (S->LeafPhoto - S->DarkResp)*((273.+DTemp)/0.53717)/(CO2 - CO2int); //eq 4 pg11
    
    // Stomatal resistance to water
    S->R_stomata = max(1e-10, 1./CndCO2 - S->RB_Water *1.3 - R_turb)/1.6;
    
    CalcPenmanMonteith(S);
    
    Dif = limit(-25.,25., S->RnetAbs - LHVAP * S->PotTran) * (R_turb + S->RB_Heat);
 
    // Second round to determine the final photosynthesis and transpiration
    InternalCO2();
    Svp_L = SatVap;
    LeafPhotoResp(S->ARAD, S->NP, &S->LeafPhoto, &S->DarkResp_n);
       
    //  Conductance for CO2
    CndCO2 = (S->LeafPhoto - S->DarkResp)*((273. + (DTemp +Dif))/0.53717)/(CO2 - CO2int); //eq 4 pg11
    
    // Adapted stomatal resistance to water
    S->R_stomata = max(1e-10, 1./CndCO2 - S->RB_Water *1.3 - R_turb)/1.6;
    
    if (Dif == 0.)
        Delta = 1.;
    else
        Delta = (Svp_L - Svp)/Dif;    //eq 6 pg 12
    
   CalcPenmanMonteith(S);
   
}

 float InstantEvap(SUSH Ev)
 {
    float SoilTemp;
    float FpeSol, FaeSol;
    float SatVap_s;    
    float PE = 0.;
    
    Dif = 0.;
    Ev->RnetAbs   = 0.;
    
    // First-round calculation to estimate soil surface temperature (Soiltemp)
    // Air-to-soil vapour pressure deficit
    SatVap     = 0.611 * exp(17.4 * DTemp / (DTemp + 239.));
    VapDeficit = max(0., SatVap - Vapour[Day][lat][lon]);
    Delta      = (4098. * SatVap)/pow((DTemp + 237.3),2); 
    
    Ev->Frac = 1.;
    Ev->R_stomata =100.;
    CalcPenmanMonteith(Ev);
            
    FpeSol = max(0., Ev->PotTran);
    FaeSol = min(FpeSol,FpeSol/(PT1+FpeSol)*WSup1);
    
    Dif = limit(-25.,25., Ev->RnetAbs - LHVAP * FaeSol) * (Ev->RT + Ev->RB_Heat);
    SoilTemp   = DTemp + Dif;

    // Second-round calculation to estimate potential soil evaporation
    SatVap_s  = 0.611 * exp(17.4 * SoilTemp / (SoilTemp + 239.));
    if (Dif == 0.)
        Delta = 1.;
    else
        Delta = (SatVap_s - SatVap)/Dif;    //eq 6 pg 12

    CalcPenmanMonteith(Ev);
          
   return max(0., PE);
    }

 void ActualAssim(SUSH S)
 {
    float ARS_Water;
    float R_turb;
    
    R_turb = RT_canopy * S->Frac;
     
    Dif = limit(-25.,25., S->RnetAbs - LHVAP * S->PotTran) * (R_turb + S->RB_Heat);
   
    // Potential photosynthesis at the new leaf temperature
    InternalCO2();
    
    // stomatal resistance to water if water stress occurs
    ARS_Water= (S->PotTran-S->ActTran)*(Delta *(S->RB_Heat +R_turb) +
            PSYCH*(S->RB_Water+R_turb))/S->ActTran/PSYCH + 
            S->PotTran/S->ActTran*S->RB_Water;

    LeafPhotoResp(S->APAR, S->NP, &S->LeafPhoto, &S->DarkResp);
    LeafPhotoResp(S->APAR, S->NP_n,&S->LeafPhoto_n,&S->DarkResp_n);

    // Actual photosynthesis under water stress condition
    S->ActualPh   = (1.6*S->R_stomata + 1.3*S->RB_Water + R_turb)/
            (1.6*ARS_Water + 1.3 * S->RB_Water + R_turb) * 
            (S->LeafPhoto - S->DarkResp) + S->DarkResp;
    S->ActualPh_n = (1.6*S->R_stomata + 1.3*S->RB_Water + R_turb)/
            (1.6*ARS_Water + 1.3 * S->RB_Water + R_turb) * 
            (S->LeafPhoto_n - S->DarkResp_n) + S->DarkResp_n;
 }


/* ----------------------------------------------------------------------------*/
/*  function DailyTotalAssimilation()                                          */ 
/*  Purpose: Calculation of the daily assimilation rate using the three point  */
/*  Gaussian integration method.                                               */
/*-----------------------------------------------------------------------------*/
float DailyTotalAssimilation()
{
    int i;
    float Hour, SinB; 
    float Rad;
    float DailyTotalAssimilation = 0.;
    float Tx, Tn;
    float GrossCO2;

    float Kl, Kwind;  // Extinction coefficient wind and nitrogen
    float Kln, Nbk;
    float PARsoil, NIRsoil;
    float FractionDiffuseRad;
    
    float BndConducH, BndConducCnp;  
    float BndConducSunlit,BndConducShaded;
    float IE, IAT, Transpiration, SoilEvap;
    float IEvap, WSup;
    float PotTran;
    
    SoilEvap  = 0.;
    Transpiration = 0.;
    
    Su = malloc(sizeof(SUSH));
    Sh = malloc(sizeof(SUSH));
    Ev = malloc(sizeof(SUSH));
    
    if (Crop->st.LAI > 0.)
    {        
        // Extinction coefficient of nitrogen  */
        Kl     = KDiff(0.2);
        Kln    = Kl * (Crop->N_st.leaves - Crop->prm.SLMIN * Crop->st.TLAI);
        Nbk    = Crop->prm.SLMIN * (1.-exp(-Kl * Crop->st.LAI));
        KNitro = 1./Crop->st.TLAI*log((Kln + Nbk)/(Kln*exp(-Kl * Crop->st.TLAI)+Nbk)); // eq41 pg35
        Kwind  =  Kl;

        // Specific leaf nitrogen and its profile in the canopy
        Sln   = Crop->N_st.leaves/Crop->st.LAI;
        Slnt  = Crop->N_st.leaves * KNitro;
        Slnnt = (Crop->N_st.leaves + 0.001*Crop->N_st.leaves)*KNitro /(1.-exp(-KNitro*Crop->st.LAI));   
        
        Fvpd = insw(Crop->C3C4, 0.195127, 0.116214); //Slope for linear effect of VPD on Ci/Ca    (kPa)-1
        
        for (i=0;i<5;i++)
        {
            Hour       = SunRise  + Daylength * XGauss[i];
            SinB       = max (0.,SinLD+CosLD*cos(2.*PI*(Hour+12.)/24.));
            Tn         = Tmin[Day][lat][lon];
            Tx         = Tmax[Day][lat][lon];
            DTemp      = Tn + (Tx-Tn) * sin(PI*(Hour + Daylength/2.-12.)/(Daylength + 3.));   
            Rad        = Radiation[Day][lat][lon] * (SinB*SolarConstant/1367.)/DSinBE;
            // Daytime course of water supply
            WSup       =  WatBal->st.RootZoneMoisture * (SinB*SolarConstant/1367.)/DSinBE;
            // Available water in the soil evaporative layer (i.e. first 5 cm)
            WSup1 = WSup * 5./Crop->st.RootDepth;
            PAR        = 0.5 * Rad;
            NIR        = PAR;
            
            if (AtmosphTransm > 0.35)
                FractionDiffuseRad = 1.47 - 1.66 * AtmosphTransm;
            if (AtmosphTransm <= 0.22 && AtmosphTransm > 0.35)
                FractionDiffuseRad = 1. - 6.4*pow((AtmosphTransm-0.22),2);
            if(AtmosphTransm < 0.22)  
                FractionDiffuseRad = 1.0;
            
            FractionDiffuseRad = max(FractionDiffuseRad, 0.15 + 0.85 * (1.-exp(-0.1/SinB)));
            
            PARDiffuse = PAR * FractionDiffuseRad;
            PARDirect  = PAR - PARDiffuse;
            NIR        = PAR;
            NIRDiffuse = PARDiffuse;
            NIRDirect  = NIR - NIRDiffuse;
            
            Kb = KBeam(SinB);
            // scattered beam radiation extinction coefficient
            KbPAR = Kb*sqrt(1.-ScatPAR);
            KbNIR = Kb*sqrt(1.-ScatNIR);
            
            KdpPAR = KDiff(ScatPAR);
            KdpNIR = KDiff(ScatNIR);

            ReflPAR_b = Reflect(ScatPAR);
            ReflNIR_b = Reflect(ScatNIR);
            
            AbsorbedLight_tot();
            
            // Photosynthetically active nitrogen for sunlit and shaded leaves
            PhotoAciveN();
            
            // Turbulence resistance for canopy (RT) and for soil (RTS) 
            Su->RT  = 0.74*pow((log((2.-0.7*Crop->st.Height) /
                    (0.1*Crop->st.Height))),2)/ (0.16 * Windspeed[Day][lat][lon]);
            Sh->RT  = Su->RT;
            Ev->RT    = 0.74*pow((log(56.)),2)/(0.16 * Windspeed[Day][lat][lon]);

            /* Canopy boundary layer conductance for sunlit and shaded leaves */
            BndConducH       = 0.01*sqrt(Windspeed[Day][lat][lon]/Crop->prm.LeafWidth);
            BndConducCnp     = (1.-exp(- 0.5 * Kwind * Crop->st.LAI))/(0.5 * Kwind ) * BndConducH;
            BndConducSunlit  = (1.-exp(-(0.5 * Kwind + Kb) * Crop->st.LAI)) /
                    (0.5 * Kwind + Kb) * BndConducH; 
            BndConducShaded  = BndConducCnp - BndConducSunlit;

            // Canopy boundary layer resistance for sunlit and shaded leaves 
            Su->RB_Heat  = 1./BndConducSunlit;       // boundary layer resistance to heat,sunlit part
            Su->RB_Water = 0.93 * Su->RB_Heat;        // boundary layer resistance to H2O, sunlit part
            Sh->RB_Heat  = 1./BndConducShaded;       // boundary layer resistance to heat,shaded part
            Sh->RB_Water = 0.93 * Sh->RB_Heat;        // boundary layer resistance to H2O, shaded part 

            // Boundary layer resistance for soil
            Ev->RB_Heat   = 172.*sqrt(0.05/max(0.1,Windspeed[Day][lat][lon]*exp(-Kwind * Crop->st.TLAI)));
            Ev->RB_Water  = 0.93*RS_Heat;
                   
            //  Fraction of sunlit leaf area 
            Su->Frac = (1./(Kb * Crop->st.LAI)) * (1 - exp(-Kb*Crop->st.LAI));
            Sh->Frac = 1 - Su->Frac;
            
            // Absorbed total radiation (PAR+NIR) by soil
            PARsoil  = 0.1;                                   //soil PAR reflection
            if (WatBal->st.Moisture - 0.5 < 0.)
                NIRsoil =  0.52-0.68 * WatBal->st.Moisture;   //soil NIR reflection
            else
                NIRsoil = 0.18;

            //  Absorbed radiation by soil
            Ev->ARAD = (1.-PARsoil)*(PARDirect * exp(-KbPAR * Crop->st.TLAI) + 
                     PARDiffuse * exp(-KdpPAR * Crop->st.TLAI)) +
                     (1.-NIRsoil)*(NIRDirect * exp(-KbNIR * Crop->st.TLAI) +
                     NIRDiffuse * exp(-KdpNIR * Crop->st.TLAI));

            InstantAssimTransp(Su);
            InstantAssimTransp(Sh);
           
            PotTran    = (Su->PotTran + Sh->PotTran);
            
            PT1 = PotTran * 5/Crop->st.RootDepth;
            // Instantaneous potential soil evaporation
            IEvap = InstantEvap(Ev);
            // Instantaneous actual soil evaporation
            IE   = min(IEvap,IEvap/(PT1+IEvap)*WSup1);
            // Instantaneous actual canopy transpiration and photosynthesis
            IAT    = min(PotTran,PT1/(PT1+IEvap)*WSup1 + WSup-WSup1);
            
            Su->ActTran = Su->PotTran/PotTran*IAT;
            Sh->ActTran = Sh->PotTran/PotTran*IAT;        
            
            // Instantaneous photosynthesis
            ActualAssim(Su);
            ActualAssim(Sh); 
            
            GrossCO2      = GrossCO2 + (Su->ActualPh + Sh->ActualPh)*WGauss[i];
            Transpiration = Transpiration + IAT*WGauss[i];
            SoilEvap      = SoilEvap + IE*WGauss[i];
            
        }  
        GrossCO2 = GrossCO2 * Daylength *3600;
        WatBal->rt.Transpiration = Transpiration * Daylength *3600;
        WatBal->rt.EvapSoil      = SoilEvap * Daylength *3600;
    }
    return(DailyTotalAssimilation*Daylength);
}


