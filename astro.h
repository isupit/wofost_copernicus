#ifndef ASTRO_H
#define ASTRO_H

float AtmosphTransm;
float AngotRadiation;
float Daylength;
float PARDaylength;
float SinLD;
float CosLD; 
float DiffRadPP;
float DSinBE;

extern float min(float a, float b);
extern float max(float a, float b);
extern float limit(float a, float b, float c);
extern float Afgen();

#endif	// ASTRO_H

