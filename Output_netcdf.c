#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "extern.h"      
#include "wofost.h"      
#include "output_netcdf.h"

/* Error handling function for NetCDF calls */
void handle_nc_error(int status) {
    if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}

/* SetupNetCDF: Creates and defines the structure of the output NetCDF file */
int SetupNetCDF(char *filename, NcFile *nc, int nlat, int nlon, int nseasons, char *tsum1_var, char *tsum2_var, char *sow_var)
{
    int dimids[3]; // To hold dimension IDs [lat, lon, year]

    /* Create the NetCDF file, overwriting if it exists */
    handle_nc_error(nc_create(filename, NC_CLOBBER | NC_NETCDF4, &nc->ncid));

    /* --- Add global attributes for traceability --- */
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum1_variable", strlen(tsum1_var), tsum1_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum2_variable", strlen(tsum2_var), tsum2_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "sow_variable", strlen(sow_var), sow_var));

    /* --- Define Dimensions --- */
    handle_nc_error(nc_def_dim(nc->ncid, "lat", nlat, &nc->lat_id));
    handle_nc_error(nc_def_dim(nc->ncid, "lon", nlon, &nc->lon_id));
    handle_nc_error(nc_def_dim(nc->ncid, "time", nseasons, &nc->time_id));    

    dimids[0] = nc->time_id;
    dimids[1] = nc->lat_id;
    dimids[2] = nc->lon_id;

    handle_nc_error(nc_def_var(nc->ncid, "Applied_N_Yearly", NC_FLOAT, 3, dimids, &nc->applied_n_yearly_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->applied_n_yearly_id, "long_name", strlen("Yearly N fertilizer application"), "Yearly N fertilizer application"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->applied_n_yearly_id, "units", strlen("kg N/ha"), "kg N/ha"));

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

    /* --- Write coordinate data --- */
    handle_nc_error(nc_put_var_double(nc->ncid, nc->lat_id, Latitude));
    handle_nc_error(nc_put_var_double(nc->ncid, nc->lon_id, Longitude));

    return NC_NOERR;
}


void WriteOutputToNetCDF(NcFile *nc)
{
    float ave, adev, sdev, var, skew, curt, lngth;
    int i;
    
    /* Start defines the [lat, lon] coordinate where we want to write */
    size_t start[2];
    start[0] = Lat;
    start[1] = Lon;
    
    /* Only write data if the simulation was successful */
    if (Crop->Seasons > 2) {

        lngth = 0;
        for (i = 1; i <= Crop->Seasons; i++) {
            lngth += Grid->length[i];
        }
        lngth /= Crop->Seasons;
        
        Moment(Grid->twso, Crop->Seasons, &ave, &adev, &sdev, &var, &skew, &curt);

        size_t start3d[3];
        size_t count3d[3];

        start3d[0] = 0;                   /* Start at the beginning of the time dimension */
        start3d[1] = Lat;                 /* Current latitude index */
        start3d[2] = Lon;                 /* Current longitude index */

        count3d[0] = Crop->Seasons;       /* Write a block of N seasons long */
        count3d[1] = 1;                   /* Write a block 1 latitude wide */
        count3d[2] = 1;                   /* Write a block 1 longitude wide */

        /* --- Write each variable to the NetCDF file --- */
        float sowing_dekad;
        int month, day;
        if (sscanf(Grid->start, "%d-%d", &month, &day) == 2) {
            int subdek = (day <= 10) ? 1 : (day <= 20) ? 2 : 3;
            sowing_dekad = (float)((month - 1) * 3 + subdek);
        } else {
            sowing_dekad = -9999.f; 
        }
        handle_nc_error(nc_put_var1_float(nc->ncid, nc->sowing_id, start, &sowing_dekad));
        
        handle_nc_error(nc_put_vara_float(nc->ncid, nc->applied_n_yearly_id, start3d, count3d, &Grid->applied_n[1]));

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

void CloseNetCDF(NcFile *nc)
{
    handle_nc_error(nc_close(nc->ncid));
}