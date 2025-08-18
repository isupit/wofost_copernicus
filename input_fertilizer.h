#ifndef INPUT_FERTILIZER_H
#define INPUT_FERTILIZER_H

#include "wofost.h" // To get the Weather struct definition

// --- Function Prototypes ---
int GetFertilizerData(Weather *meteo, char *filename, char *varname);
void CleanFertilizerData(Weather* meteo);

#endif