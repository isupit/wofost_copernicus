#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "extern.h"      // Assuming this contains Latitude, Longitude etc.
#include "wofost.h"      // Assuming this contains Crop, Grid, Moment() etc.
#include "output_netcdf.h"

/*
 * A simple error handling function for NetCDF calls.
 */
void handle_nc_error(int status) {
    if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}

/*
 * SETUPNETCDF: Creates and defines the structure of the output NetCDF file.
 */
int SetupNetCDF(char *filename, NcFile *nc, int nlat, int nlon, char *tsum1_var, char *tsum2_var, char *sow_var)
{
    int dimids[2]; // To hold dimension IDs [lat, lon]

    /* Create the NetCDF file, overwriting if it exists */
    handle_nc_error(nc_create(filename, NC_CLOBBER | NC_NETCDF4, &nc->ncid));

    /* --- Add global attributes for traceability --- */
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum1_variable", strlen(tsum1_var), tsum1_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum2_variable", strlen(tsum2_var), tsum2_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "sow_variable", strlen(sow_var), sow_var));

    /* --- Define Dimensions --- */
    handle_nc_error(nc_def_dim(nc->ncid, "lat", nlat, &nc->lat_id));
    handle_nc_error(nc_def_dim(nc->ncid, "lon", nlon, &nc->lon_id));
    
    dimids[0] = nc->lat_id;
    dimids[1] = nc->lon_id;

    /* --- Define Coordinate Variables (lat, lon) --- */
    handle_nc_error(nc_def_var(nc->ncid, "lat", NC_DOUBLE, 1, &nc->lat_id, &nc->lat_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->lat_id, "units", strlen("degrees_north"), "degrees_north"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->lat_id, "long_name", strlen("latitude"), "latitude"));

    handle_nc_error(nc_def_var(nc->ncid, "lon", NC_DOUBLE, 1, &nc->lon_id, &nc->lon_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->lon_id, "units", strlen("degrees_east"), "degrees_east"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->lon_id, "long_name", strlen("longitude"), "longitude"));

    /* --- Define Data Variables (all are 2D with dimensions lat, lon) --- */
    handle_nc_error(nc_def_var(nc->ncid, "SowingDate", NC_FLOAT, 2, dimids, &nc->sowing_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->sowing_id, "long_name", strlen("Sowing date"), "Sowing date"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->sowing_id, "units", strlen("dekad"), "dekad"));
    handle_nc_error(nc_put_att_float(nc->ncid, nc->sowing_id, "_FillValue", NC_FLOAT, 1, &(float){-9999.f}));

    handle_nc_error(nc_def_var(nc->ncid, "Length", NC_FLOAT, 2, dimids, &nc->length_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->length_id, "long_name", strlen("Average length of growing season"), "Average length of growing season"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->length_id, "units", strlen("days"), "days"));

    handle_nc_error(nc_def_var(nc->ncid, "TSM1", NC_FLOAT, 2, dimids, &nc->tsm1_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->tsm1_id, "long_name", strlen("Temperature sum emergence to anthesis"), "Temperature sum emergence to anthesis"));
    
    handle_nc_error(nc_def_var(nc->ncid, "TSM2", NC_FLOAT, 2, dimids, &nc->tsm2_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->tsm2_id, "long_name", strlen("Temperature sum anthesis to maturity"), "Temperature sum anthesis to maturity"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Average", NC_FLOAT, 2, dimids, &nc->avg_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->avg_id, "units", strlen("kg/ha"), "kg/ha"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Avg_Deviation", NC_FLOAT, 2, dimids, &nc->adev_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->adev_id, "units", strlen("kg/ha"), "kg/ha"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Std_Deviation", NC_FLOAT, 2, dimids, &nc->sdev_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->sdev_id, "units", strlen("kg/ha"), "kg/ha"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Variance", NC_FLOAT, 2, dimids, &nc->var_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->var_id, "units", strlen("kg2/ha2"), "kg2/ha2"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Skewness", NC_FLOAT, 2, dimids, &nc->skew_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->skew_id, "units", strlen("dimensionless"), "dimensionless"));

    handle_nc_error(nc_def_var(nc->ncid, "Yield_Kurtosis", NC_FLOAT, 2, dimids, &nc->curt_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->curt_id, "units", strlen("dimensionless"), "dimensionless"));

    handle_nc_error(nc_def_var(nc->ncid, "Simulated_Seasons", NC_INT, 2, dimids, &nc->seasons_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->seasons_id, "long_name", strlen("Number of simulated seasons"), "Number of simulated seasons"));

    /* End define mode */
    handle_nc_error(nc_enddef(nc->ncid));

    /* --- Write coordinate data (since they don't change) --- */
    handle_nc_error(nc_put_var_double(nc->ncid, nc->lat_id, Latitude));
    handle_nc_error(nc_put_var_double(nc->ncid, nc->lon_id, Longitude));

    return NC_NOERR;
}

/*
 * WRITEOUTPUTTONETCDF: Calculates values and writes them for the current grid cell.
 */
void WriteOutputToNetCDF(NcFile *nc)
{
    float ave, adev, sdev, var, skew, curt, lngth;
    int i;
    
    /* start defines the [lat, lon] coordinate where we want to write */
    size_t start[2];
    start[0] = Lat;
    start[1] = Lon;
    
    /* Only write data if the simulation was successful */
    if (Crop->Seasons > 2) {
        /* --- Perform Calculations --- */
        
        lngth = 0;
        for (i = 1; i <= Crop->Seasons; i++) {
            lngth += Grid->length[i];
        }
        lngth /= Crop->Seasons;
        
        Moment(Grid->twso, Crop->Seasons, &ave, &adev, &sdev, &var, &skew, &curt);

        /* --- Write each variable to its place in the NetCDF file --- */
        
        // Convert "MM-DD" from Grid->start back to dekad (float)
        float sowing_dekad;
        int month, day;
        if (sscanf(Grid->start, "%d-%d", &month, &day) == 2) {
            // Approximate dekad: (month-1)*3 + ceil(day/10)
            int subdek = (day <= 10) ? 1 : (day <= 20) ? 2 : 3;
            sowing_dekad = (float)((month - 1) * 3 + subdek);
        } else {
            sowing_dekad = -9999.f; // Use _FillValue for invalid format
        }
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->sowing_id, start, &sowing_dekad));
        
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->length_id, start, &lngth));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->tsm1_id, start, &Crop->prm.TempSum1));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->tsm2_id, start, &Crop->prm.TempSum2));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->avg_id, start, &ave));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->adev_id, start, &adev));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->sdev_id, start, &sdev));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->var_id, start, &var));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->skew_id, start, &skew));
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->curt_id, start, &curt));
        handle_nc_error(nc_put_var1_int(nc->ncid, nc->seasons_id, start, &Crop->Seasons));
    }
}

/*
 * CLOSENETCDF: Closes the NetCDF file.
 */
void CloseNetCDF(NcFile *nc)
{
    handle_nc_error(nc_close(nc->ncid));
}