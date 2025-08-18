#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <netcdf.h>
#include "wofost.h"
#include "extern.h"
#include "input_fertilizer.h"

// Simple error handler for NetCDF calls
void handle_nc_fert_error(int status) {
    if (status != NC_NOERR) {
        fprintf(stderr, "NetCDF Error: %s\n", nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}


int GetFertilizerData(Weather *meteo, char *filename, char *varname) {
    int retval, ncid, varid;
    int lat_dimid, lon_dimid, time_dimid;
    size_t lat_len, lon_len;
    
    // Assumes time units are "days since YYYY-01-01 00:00:00" or similar
    char time_units[MAX_STRING];
    int start_year;

    handle_nc_fert_error(nc_open(filename, NC_NOWRITE, &ncid));

    // --- Get dimensions and check for consistency ---
    handle_nc_fert_error(nc_inq_dimid(ncid, "lat", &lat_dimid));
    handle_nc_fert_error(nc_inq_dimid(ncid, "lon", &lon_dimid));
    handle_nc_fert_error(nc_inq_dimid(ncid, "year", &time_dimid));

    handle_nc_fert_error(nc_inq_dimlen(ncid, lat_dimid, &lat_len));
    handle_nc_fert_error(nc_inq_dimlen(ncid, lon_dimid, &lon_len));
    handle_nc_fert_error(nc_inq_dimlen(ncid, time_dimid, &meteo->n_fert_time_len));

    if (lat_len != meteo->nlat || lon_len != meteo->nlon) {
        fprintf(stderr, "Fertilizer grid dimensions (%zu, %zu) do not match meteo grid (%zu, %zu).\n",
                lat_len, lon_len, meteo->nlat, meteo->nlon);
        exit(1);
    }

    meteo->n_fert_start_year = 1961;
    printf("Set fertilizer data start year to %d (hardcoded).\n", meteo->n_fert_start_year);

    // --- Allocate memory for the fertilizer grid [time][lat][lon] ---
    meteo->n_fertilizer_grid = malloc(meteo->n_fert_time_len * sizeof(*meteo->n_fertilizer_grid));
    for (size_t t = 0; t < meteo->n_fert_time_len; t++) {
        meteo->n_fertilizer_grid[t] = malloc(lat_len * sizeof(*meteo->n_fertilizer_grid[t]));
        for (size_t j = 0; j < lat_len; j++) {
            meteo->n_fertilizer_grid[t][j] = malloc(lon_len * sizeof(*meteo->n_fertilizer_grid[t][j]));
        }
    }

    // --- Read the data ---
    handle_nc_fert_error(nc_inq_varid(ncid, varname, &varid));

    printf("Reading N fertilizer variable '%s' from %s\n", varname, filename);
    
    // Read all data at once for efficiency
    float *temp_data = malloc(meteo->n_fert_time_len * lat_len * lon_len * sizeof(float));
    handle_nc_fert_error(nc_get_var_float(ncid, varid, temp_data));

    // Distribute the flat array into the 3D pointer structure
    for (size_t t = 0; t < meteo->n_fert_time_len; t++) {
        for (size_t j = 0; j < lat_len; j++) {
            for (size_t k = 0; k < lon_len; k++) {
                meteo->n_fertilizer_grid[t][j][k] = temp_data[t * (lat_len * lon_len) + j * lon_len + k];
            }
        }
    }
    
    free(temp_data);
    handle_nc_fert_error(nc_close(ncid));

    return 1;
}

// Function to free the allocated memory
void CleanFertilizerData(Weather* meteo) {
    if (meteo->n_fertilizer_grid) {
        for (size_t t = 0; t < meteo->n_fert_time_len; t++) {
            for (size_t j = 0; j < meteo->nlat; j++) {
                free(meteo->n_fertilizer_grid[t][j]);
            }
            free(meteo->n_fertilizer_grid[t]);
        }
        free(meteo->n_fertilizer_grid);
        meteo->n_fertilizer_grid = NULL;
    }
}