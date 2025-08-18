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

int SetupNetCDF(char *filename, NcFile *nc, int nlat, int nlon, int nseasons,
                char *tsum1_var, char *tsum2_var, char *sow_var)
{
    /* Create file */
    handle_nc_error(nc_create(filename, NC_CLOBBER | NC_NETCDF4, &nc->ncid));

    /* Global attrs */
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum1_variable", strlen(tsum1_var), tsum1_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "tsum2_variable", strlen(tsum2_var), tsum2_var));
    handle_nc_error(nc_put_att_text(nc->ncid, NC_GLOBAL, "sow_variable",   strlen(sow_var),   sow_var));

    /* --- Define dimensions --- */
    int dim_lat_id, dim_lon_id, dim_time_id;
    handle_nc_error(nc_def_dim(nc->ncid, "lat",  nlat,     &dim_lat_id));
    handle_nc_error(nc_def_dim(nc->ncid, "lon",  nlon,     &dim_lon_id));
    handle_nc_error(nc_def_dim(nc->ncid, "time", nseasons, &dim_time_id));

    /* Keep separate dim arrays */
    int dims_latlon[2]      = { dim_lat_id, dim_lon_id };
    int dims_timelatlon[3]  = { dim_time_id, dim_lat_id, dim_lon_id };

    /* --- Coordinate vars --- */
    int lat_varid, lon_varid;
    handle_nc_error(nc_def_var(nc->ncid, "lat", NC_DOUBLE, 1, &dim_lat_id, &lat_varid));
    handle_nc_error(nc_put_att_text(nc->ncid, lat_varid, "units", strlen("degrees_north"), "degrees_north"));
    handle_nc_error(nc_put_att_text(nc->ncid, lat_varid, "long_name", strlen("latitude"), "latitude"));

    handle_nc_error(nc_def_var(nc->ncid, "lon", NC_DOUBLE, 1, &dim_lon_id, &lon_varid));
    handle_nc_error(nc_put_att_text(nc->ncid, lon_varid, "units", strlen("degrees_east"), "degrees_east"));
    handle_nc_error(nc_put_att_text(nc->ncid, lon_varid, "long_name", strlen("longitude"), "longitude"));

    /* store varids for later writes */
    nc->lat_id  = lat_varid;
    nc->lon_id  = lon_varid;
    nc->time_id = dim_time_id; /* only used as a dim id later */

    /* --- Data vars --- */
    /* 3-D: time,lat,lon */
    handle_nc_error(nc_def_var(nc->ncid, "Applied_N_Yearly", NC_FLOAT, 3, dims_timelatlon, &nc->applied_n_yearly_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->applied_n_yearly_id, "long_name",
                                    strlen("Yearly N fertilizer application"), "Yearly N fertilizer application"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->applied_n_yearly_id, "units", strlen("kg N/ha"), "kg N/ha"));

    /* 2-D: lat,lon  (IMPORTANT: use dims_latlon, not the 3-D array) */
    handle_nc_error(nc_def_var(nc->ncid, "SowingDate", NC_FLOAT, 2, dims_latlon, &nc->sowing_id));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->sowing_id, "long_name", strlen("Sowing date"), "Sowing date"));
    handle_nc_error(nc_put_att_text(nc->ncid, nc->sowing_id, "units", strlen("dekad"), "dekad"));
    handle_nc_error(nc_put_att_float(nc->ncid, nc->sowing_id, "_FillValue", NC_FLOAT, 1, &(float){-9999.f}));

    handle_nc_error(nc_def_var(nc->ncid, "Length",          NC_FLOAT, 2, dims_latlon, &nc->length_id));
    handle_nc_error(nc_def_var(nc->ncid, "TSM1",            NC_FLOAT, 2, dims_latlon, &nc->tsm1_id));
    handle_nc_error(nc_def_var(nc->ncid, "TSM2",            NC_FLOAT, 2, dims_latlon, &nc->tsm2_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Average",   NC_FLOAT, 2, dims_latlon, &nc->avg_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Avg_Deviation", NC_FLOAT, 2, dims_latlon, &nc->adev_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Std_Deviation", NC_FLOAT, 2, dims_latlon, &nc->sdev_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Variance",      NC_FLOAT, 2, dims_latlon, &nc->var_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Skewness",      NC_FLOAT, 2, dims_latlon, &nc->skew_id));
    handle_nc_error(nc_def_var(nc->ncid, "Yield_Kurtosis",      NC_FLOAT, 2, dims_latlon, &nc->curt_id));
    handle_nc_error(nc_def_var(nc->ncid, "Simulated_Seasons",   NC_INT,   2, dims_latlon, &nc->seasons_id));

    /* End define mode and write coords */
    handle_nc_error(nc_enddef(nc->ncid));
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