#ifndef ASSIM_H
#define ASSIM_H

#define BOLTZM  5.668e-8       // Stefan-Boltzmann constant(J/m2/s/K4)
#define LHVAP   2.45e6           // Latent heat of water vaporization(J/kg)
#define VHCA    1200.           // Volumetric heat capacity (J/m3/oC)
#define PSYCH   0.067           // Psychrometric constant (kPa/oC)


float ReflPAR_b,ReflNIR_b;
float ReflPAR_d,ReflNIR_d;

float PARDirect, NIRDirect;
float PARDiffuse,NIRDiffuse;

float KbPAR, KbNIR;
float KdpPAR,KdpNIR;

float PAR, NIR;

float ASoil;

float ScatP;
float ScatN;
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
float Dif;
float KdiffPAR, KdiffNIR;
float AtmTransm;
float SinB;
float Fvpd;
float SatVap, VapDeficit, Delta;
float PT1;
float WSup1;
#endif	// ASSIM_H

