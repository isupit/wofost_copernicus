
#ifndef PENMAN_H
#define PENMAN_H

typedef struct ETP {
    float E0;
    float ES0;
    float ET0;
} Etp;

Etp Penman;

extern float min(float a, float b);
extern float max(float a, float b);
extern float limit(float a, float b, float c);

typedef struct EVP {
    float MaxEvapWater;
    float MaxEvapSoil;
    float MaxTranspiration;
} EVP;
EVP Evtra;

float RB_Heat_su;
float RB_Water_su;
float RB_Heat_sh; 
float RB_Water_sh;

float RT_canopy;    // Turbulence resistance canopy
float RT_soil;      // Turbulence resistance soil

float RS_Heat;     // Boundary resistance soil heat
float RS_Water;    // Boundary resistance soil water

float R_stomata;   //Stomatal resistance

float FracSunlit;
float FracShaded;

float APARDiffuse, APARTotal, APARDirect;
float ANIRDiffuse, ANIRTotal, ANIRDirect;
float APARShadedLeaves, APARSunlitLeaves;
float ANIRShadedLeaves, ANIRSunlitLeaves;
float AShadedLeaves, ASunlitLeaves;
float ASoil;
float AssimShaded, AssimSunlit;


#endif	// PENMAN_H

