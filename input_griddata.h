#ifndef INPUT_GRIDDATA_H
#define INPUT_GRIDDATA_H

#include "wofost.h"

/* --- Function Prototypes --- */

/* Reads all gridded input data (TSum, Sowing Date, etc.) from a NetCDF file */
void GetGridData(Weather *meteo, char *grid_nc_file, char *tsum1_var, char *tsum2_var, char *sow_var);

/* Frees the memory allocated for all gridded input data */ 
void CleanGridData(Weather *meteo); 

#endif /* INPUT_GRIDDATA_H */ 