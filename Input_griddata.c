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

void SetDefaultTSUM(Weather *meteo)
{
    size_t i, j;
    
    // Get default TSUM values from the first crop parameter (assuming all crops have same defaults)
    SimUnit *tempGrid = Grid;
    if (tempGrid == NULL || tempGrid->crp == NULL) {
        fprintf(stderr, "Error: No crop parameters available for default TSUM\n");
        exit(1);
    }
    float default_tsum1 = tempGrid->crp->prm.TempSum1;
    float default_tsum2 = tempGrid->crp->prm.TempSum2;
    
    printf("Setting default TSUM1=%.1f, TSUM2=%.1f for all %zux%zu grid cells\n", 
           default_tsum1, default_tsum2, meteo->nlat, meteo->nlon);
    
    // Allocate memory for 2D crop arrays in the Meteo struct if not already allocated
    if (meteo->tsum1_grid == NULL) {
        meteo->tsum1_grid = malloc(meteo->nlat * sizeof(float *));
        meteo->tsum2_grid = malloc(meteo->nlat * sizeof(float *));
        meteo->sowing_date_grid = malloc(meteo->nlat * sizeof(float *));  // Allocate even for defaults
        
        if (meteo->tsum1_grid == NULL || meteo->tsum2_grid == NULL || meteo->sowing_date_grid == NULL) {
            fprintf(stderr, "Error: Failed to allocate memory for default TSUM grids\n");
            exit(1);
        }
        
        for (i = 0; i < meteo->nlat; i++) {
            meteo->tsum1_grid[i] = malloc(meteo->nlon * sizeof(float));
            meteo->tsum2_grid[i] = malloc(meteo->nlon * sizeof(float));
            meteo->sowing_date_grid[i] = malloc(meteo->nlon * sizeof(float));
            
            if (meteo->tsum1_grid[i] == NULL || meteo->tsum2_grid[i] == NULL || meteo->sowing_date_grid[i] == NULL) {
                fprintf(stderr, "Error: Failed to allocate memory for row %zu in default TSUM grids\n", i);
                exit(1);
            }
        }
    }
    
    // Set default values for all grid cells
    for (i = 0; i < meteo->nlat; i++) {
        for (j = 0; j < meteo->nlon; j++) {
            meteo->tsum1_grid[i][j] = default_tsum1;
            meteo->tsum2_grid[i][j] = default_tsum2;
            meteo->sowing_date_grid[i][j] = 10.0f;  // Default sowing date: April 1st (dekad 10)
        }
    }
}

void LoadSowingDateOnly(Weather *meteo, char *grid_nc_file, char *sow_var)
{
    int ncid, sow_date_id; 
    size_t i, j;
    double *temp_sow_dates; 
    
    printf("Loading only sowing date data from %s (%s variable)...\n", grid_nc_file, sow_var);
    
    /* Open the NetCDF file */
    handle_grid_nc_error(nc_open(grid_nc_file, NC_NOWRITE, &ncid));
    
    /* Get variable ID for sowing date */
    handle_grid_nc_error(nc_inq_varid(ncid, sow_var, &sow_date_id));
    
    /* Allocate memory for sowing date grid if not already allocated */
    if (meteo->sowing_date_grid == NULL) {
        meteo->sowing_date_grid = malloc(meteo->nlat * sizeof(float *));
        for (i = 0; i < meteo->nlat; i++) {
            meteo->sowing_date_grid[i] = malloc(meteo->nlon * sizeof(float));
        }
    }
    
    /* Allocate memory for temporary flat array */
    temp_sow_dates = malloc(meteo->nlat * meteo->nlon * sizeof(double));
    
    /* Read the sowing date data */
    handle_grid_nc_error(nc_get_var_double(ncid, sow_date_id, temp_sow_dates)); 
    
    /* Copy data from flat temp array to 2D array */
    for (i = 0; i < meteo->nlat; i++) {
        for (j = 0; j < meteo->nlon; j++) {
            meteo->sowing_date_grid[i][j] = (float)temp_sow_dates[i * meteo->nlon + j];
        }
    }
    
    /* Clean up temporary array */
    free(temp_sow_dates);
    
    handle_grid_nc_error(nc_close(ncid));
    printf("Sowing date data loaded successfully.\n");
}

/* Apply TSUM offsets to all grid cells */
void ApplyTSUMOffsets(Weather *meteo)
{
    size_t i, j;
    
    // FIXED: Declare offset variables as extern since they're defined in Wofost.c
    extern float tsum1_offset;
    extern float tsum2_offset;
    
    for (i = 0; i < meteo->nlat; i++) {
        for (j = 0; j < meteo->nlon; j++) {
            meteo->tsum1_grid[i][j] += tsum1_offset;
            meteo->tsum2_grid[i][j] += tsum2_offset;
        }
    }
    
    // FIXED: Use %zu for size_t
    printf("Applied offsets to all %zux%zu grid cells successfully.\n", meteo->nlat, meteo->nlon);
}
