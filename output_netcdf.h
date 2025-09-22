#ifndef OUTPUT_NETCDF_H
#define OUTPUT_NETCDF_H

#include <stdio.h>
#include <netcdf.h>

/* This struct will hold the NetCDF IDs for the file and all variables. */
typedef struct {
    int ncid;
    int lat_id, lon_id, time_id;
    int sowing_id, length_id, tsm1_id, tsm2_id;
    int avg_id, adev_id, sdev_id, var_id;
    int skew_id, curt_id, seasons_id;
    int applied_n_yearly_id; 
    int cold_days_yearly_id;   
    int n_storage_id, p_storage_id, k_storage_id;
} NcFile;

/* --- Function Prototypes --- */
/* Initializes the NetCDF file, defines dimensions and variables. */
int SetupNetCDF(char *filename, NcFile *nc, int nlat, int nlon, int nseasons, char *tsum1_var, char *tsum2_var, char *sow_var);
/* Writes the output for a single grid cell to the NetCDF file. */
void WriteOutputToNetCDF(NcFile *nc);

/* Closes the NetCDF file. */
void CloseNetCDF(NcFile *nc);

#endif /* OUTPUT_NETCDF_H */
