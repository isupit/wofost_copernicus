#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "wofost.h"
#include "extern.h"
#include "input_griddata.h" 

/* A simple error handling function for NetCDF calls. */
void handle_grid_nc_error(int status) { // <<< MODIFIED: Renamed for clarity
    if (status != NC_NOERR) {
        fprintf(stderr, "NetCDF Error (Grid Input): %s\n", nc_strerror(status));
        exit(EXIT_FAILURE);
    }
}

/*
 * GETGRIDDATA: Reads TSM and Sowing Date data from the specified NetCDF file
 * and stores them in the Meteo structure.
 */
void GetGridData(Weather *meteo, char *grid_nc_file, char *tsum1_var, char *tsum2_var, char *sow_var)
{
    int ncid, tsum1_id, tsum2_id, sow_date_id; 
    size_t lat_len, lon_len;
    double *temp_tsum1, *temp_tsum2;
    double *temp_sow_dates; 
    size_t i, j;

    printf("Reading gridded data from %s...\n", grid_nc_file);

    /* Open the NetCDF file */
    handle_grid_nc_error(nc_open(grid_nc_file, NC_NOWRITE, &ncid));

    /* Get and verify dimensions (this part remains the same) */
    // ... (code for nc_inq_dimid, nc_inq_dimlen, and dimension check) ...

    /* Get variable IDs */
    handle_grid_nc_error(nc_inq_varid(ncid, tsum1_var, &tsum1_id));
    handle_grid_nc_error(nc_inq_varid(ncid, tsum2_var, &tsum2_id));
    handle_grid_nc_error(nc_inq_varid(ncid, sow_var, &sow_date_id));

    /* Allocate memory for 2D crop arrays in the Meteo struct */
    meteo->tsum1_grid = malloc(meteo->nlat * sizeof(float *));
    meteo->tsum2_grid = malloc(meteo->nlat * sizeof(float *));
    meteo->sowing_date_grid = malloc(meteo->nlat * sizeof(float *));

    for (i = 0; i < meteo->nlat; i++) {
        meteo->tsum1_grid[i] = malloc(meteo->nlon * sizeof(float));
        meteo->tsum2_grid[i] = malloc(meteo->nlon * sizeof(float));
        meteo->sowing_date_grid[i] = malloc(meteo->nlon * sizeof(float));
    }
    
    /* Allocate memory for temporary flat arrays to read into */
    temp_tsum1 = malloc(meteo->nlat * meteo->nlon * sizeof(double));
    temp_tsum2 = malloc(meteo->nlat * meteo->nlon * sizeof(double));
    // +++ MODIFIED: Use double* for sowing dates (dekads are float)
    temp_sow_dates = malloc(meteo->nlat * meteo->nlon * sizeof(double));

    /* Read the data */
    handle_grid_nc_error(nc_get_var_double(ncid, tsum1_id, temp_tsum1));
    handle_grid_nc_error(nc_get_var_double(ncid, tsum2_id, temp_tsum2));
    handle_grid_nc_error(nc_get_var_double(ncid, sow_date_id, temp_sow_dates)); 

    /* Copy data from flat temp arrays to 2D arrays */
    for (i = 0; i < meteo->nlat; i++) {
        for (j = 0; j < meteo->nlon; j++) {
            meteo->tsum1_grid[i][j] = (float)temp_tsum1[i * meteo->nlon + j];
            meteo->tsum2_grid[i][j] = (float)temp_tsum2[i * meteo->nlon + j];
            meteo->sowing_date_grid[i][j] = (float)temp_sow_dates[i * meteo->nlon + j];
        }
    }
    
    /* Clean up temporary arrays */
    free(temp_tsum1);
    free(temp_tsum2);
    free(temp_sow_dates);
    
    handle_grid_nc_error(nc_close(ncid));
    printf("Spatially-variable grid data loaded successfully.\n");
}


void CleanGridData(Weather *meteo) 
{
    size_t i;
    // Tsum cleanup
    if (meteo->tsum1_grid != NULL) { 
        for (i = 0; i < meteo->nlat; i++) {
            free(meteo->tsum1_grid[i]);
        }
        free(meteo->tsum1_grid);
        meteo->tsum1_grid = NULL;
    }
    if (meteo->tsum2_grid != NULL) { 
        for (i = 0; i < meteo->nlat; i++) {
            free(meteo->tsum2_grid[i]);
        }
        free(meteo->tsum2_grid);
        meteo->tsum2_grid = NULL;
    }

    // Free the sowing date grid 
    if (meteo->sowing_date_grid != NULL) {
        for (i = 0; i < meteo->nlat; i++) {
            free(meteo->sowing_date_grid[i]);
        }
        free(meteo->sowing_date_grid);
        meteo->sowing_date_grid = NULL;
    }
}