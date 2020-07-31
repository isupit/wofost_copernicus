#ifndef ASSIM_H
#define ASSIM_H

#define BOLTZM  5.668 E-8       // Stefan-Boltzmann constant(J/m2/s/K4)
#define LHVAP   2.45e6           // Latent heat of water vaporization(J/kg)
#define VHCA    1200.           // Volumetric heat capacity (J/m3/oC)
#define PSYCH   0.067           // Psychrometric constant (kPa/oC)

float DayTemp;
float VapourDef;
float SatVap;
float Slope;
float Dif;
float Kb, Kw;
float KdiffPAR, KdiffNIR;
float PARDiffuse, PARDirect;
float NIRDiffuse, NIRDirect;
float AtmTransm;
float FracSunlit, FracShaded;
float Rtc, Rts;
float SnlTotAbsorbed, ShTotAbsorbed;
float SinB;
float Fvpd;
float SatVap, VapDeficit, Delta;
float LeafPhoto, PotTran;





#endif	// ASSIM_H

