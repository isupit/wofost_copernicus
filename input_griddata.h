#ifndef INPUT_GRIDDATA_H
#define INPUT_GRIDDATA_H

#include "wofost.h"

/* --- Function Prototypes --- */

/* * Reads all gridded input data (TSM, Sowing Date, etc.) from a NetCDF file. 
 */ // <<< MODIFIED: Updated comment
void GetGridData(Weather *meteo, char *grid_nc_file, char *tsum1_var, char *tsum2_var, char *sow_var);

/* * Frees the memory allocated for all gridded input data.
 */ // <<< MODIFIED: Updated comment
void CleanGridData(Weather *meteo); // <<< MODIFIED: Renamed function

#endif /* INPUT_GRIDDATA_H */ // <<< MODIFIED: Renamed header guard