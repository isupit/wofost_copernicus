#ifndef ASSIM_H
#define ASSIM_H

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

#define BOLTZM  5.668e-8        // Stefan-Boltzmann constant(J/m2/s/K4)
#define LHVAP   2.45e6          // Latent heat of water vaporization(J/kg)
#define VHCA    1200.           // Volumetric heat capacity (J/m3/oC)
#define PSYCH   0.067           // Psychrometric constant (kPa/oC)

float ReflPAR_b,ReflNIR_b;
float PARDirect, NIRDirect;
float PARDiffuse,NIRDiffuse;

float KbPAR, KbNIR;
float KdpPAR,KdpNIR;

float ASoil;
float Kb, Kw;
float Kbp, Kdp;
        
float KNitro, Sln, Slnt, Slnnt;
float RT_canopy, RT_soil;

typedef struct SUSH {
    float Frac;
    float NP, NP_n;
    float APAR, ANIR;
    float ARAD;
    float RT;
    float RB_Heat;
    float RB_Water;
    float R_stomata;
    float PotTran;
    float ActTran;
    float LeafPhoto;
    float LeafPhoto_n;
    float ActualPh;
    float ActualPh_n;
    float DarkResp;
    float DarkResp_n;
    float RnetAbs;
} *SUSH;
SUSH Su, Sh, Ev;

float RS_Heat, RS_Water;

float DTemp;
float CO2int;

float Slope;
float Dif;   // Temperature difference plant/soil and air temperature.
float AtmTransm;
float SinB;
float Fvpd;
float SatVap, VapDeficit, Delta;
float PT1;
float WSup1;
#endif	// ASSIM_H

