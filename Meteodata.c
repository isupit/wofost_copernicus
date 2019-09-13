#include <stdio.h>#include <string.h>#include <stdlib.h>#include <math.h>#include <netcdf.h>#include "wofost.h"/* Handle errors by printing an error message and exiting with a * non-zero status. */#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}/* set decimals */#define roundz(x,d) ((floor(((x)*pow(10,d))+.5))/pow(10,d))int GetMeteoData(char* path, char* inputfile){    int i;    size_t lt, ln;    int ncid, varid;    int lat_varid, lon_varid, time_varid;           /* netcdf input files */    char prec_name[MAX_STRING];    char hurs_name[MAX_STRING];    char rsds_name[MAX_STRING];    char tmax_name[MAX_STRING];    char tmin_name[MAX_STRING];    char wind_name[MAX_STRING];        /* Error handling. */    int retval;        //float temp;    float Svap_Tmin, Svap_Tmax, Svap;            AngstA = 0.249;    AngstB = 0.515;    Altitude = 1200;        /* Construct the rainfall file name */    memset(prec_name,'\0',MAX_STRING);    if ((strlen(inputfile) + strlen(path) + 10) > MAX_STRING) exit(0);    strncpy(prec_name, path,strlen(path));    strcat(prec_name,"pr");    strncat(prec_name, inputfile, strlen(inputfile));        /* Construct the relative humidity file file name */    memset(hurs_name,'\0',MAX_STRING);    strncpy(hurs_name, path,strlen(path));    strcat(hurs_name,"hurs");    strncat(hurs_name, inputfile, strlen(inputfile));    /* Construct the radiation file name */    memset(rsds_name,'\0',MAhttp://www.fao.org/3/X0490E/x0490e06.htm#fao%20penman%20monteith%20equationX_STRING);    strncpy(rsds_name, path,strlen(path));    strcat(rsds_name,"rsds");    strncat(rsds_name, inputfile, strlen(inputfile));        /* Construct the tmax file name */    memset(tmax_name,'\0',MAX_STRING);    strncpy(tmax_name, path,strlen(path));    strcat(tmax_name,"tasmax");    strncat(tmax_name, inputfile, strlen(inputfile));         /* Construct the tmin file name */    memset(tmin_name,'\0',MAX_STRING);    strncpy(tmin_name, path,strlen(path));    strcat(tmin_name,"tasmin");    strncat(tmin_name, inputfile, strlen(inputfile));          /* Construct the wind file name */    memset(wind_name,'\0',MAX_STRING);    strncpy(wind_name, path,strlen(path));    strcat(wind_name,"wind");    strncat(wind_name, inputfile, strlen(inputfile));         /* Open the precipitation file. */    if ((retval = nc_open(prec_name, NC_NOWRITE, &ncid)))       ERR(retval);        /* Get the varids of rain, time, latitude and longitude variables */    if ((retval = nc_inq_varid(ncid, "prAdjust", &varid)))        ERR(retval);    if((retval  = nc_inq_dimid(ncid, "time", &time_varid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "lat", &lat_varid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "lon", &lon_varid)))        ERR(retval);        /* Get the length of the time, lats and lons */    if ((retval = nc_inq_dimlen(ncid, time_varid, &time_length)))       ERR(retval);        if ((retval = nc_inq_dimlen(ncid, lat_varid, &lat_length)))       ERR(retval);     if ((retval = nc_inq_dimlen(ncid, lon_varid, &lon_length)))       ERR(retval);        /* Get the values of the rain, lat, lon variables */    /* Note that the values of the time values are read ! */    if ((retval = nc_get_var_float(ncid, lon_varid, &lons[0])))       ERR(retval);    if ((retval = nc_get_var_float(ncid, lat_varid, &lats[0])))       ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Rain[0][0][0])))       ERR(retval);    if ((retval = nc_close(ncid)))      ERR(retval);        /* Open the radiation file. */    if ((retval = nc_open(rsds_name, NC_NOWRITE, &ncid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "rsdsAdjust", &varid)))        ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Radiation[0][0][0])))        ERR(retval);    if ((retval = nc_close(ncid)))        ERR(retval);            /* Open the relative humidity file. */    if ((retval = nc_open(hurs_name, NC_NOWRITE, &ncid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "hurs", &varid)))        ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Vapour[0][0][0])))        ERR(retval);    if ((retval = nc_close(ncid)))        ERR(retval);        /* Open the maximum temperature file. */    if ((retval = nc_open(tmax_name, NC_NOWRITE, &ncid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "tasmaxAdjust", &varid)))        ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Tmax[0][0][0])))        ERR(retval);    if ((retval = nc_close(ncid)))        ERR(retval);        /* Open the minimum temperature file. */    if ((retval = nc_open(tmin_name, NC_NOWRITE, &ncid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "tasminAdjust", &varid)))        ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Tmin[0][0][0])))        ERR(retval);    if ((retval = nc_close(ncid)))        ERR(retval);        /* Open the wind file. */    if ((retval = nc_open(wind_name, NC_NOWRITE, &ncid)))        ERR(retval);    if ((retval = nc_inq_varid(ncid, "windAdjust", &varid)))        ERR(retval);    if ((retval = nc_get_var_float(ncid, varid, &Windspeed[0][0][0])))        ERR(retval);    if ((retval = nc_close(ncid)))        ERR(retval);            for (lt = 0; lt < lat_length; lt++)    {          for (ln = 0; ln < lon_length; ln++)        {               for (i = time_length - 1; i >= 0; i-- )            {                Rain[i + 1][lt][ln] = roundz(8640 * Rain[i][lt][ln],2);   // from kg m-1 sec to cm day-1                Tmax[i + 1][lt][ln] = roundz(Tmax[i][lt][ln] - 273.15,1); // from Kelvin to Celsius                Tmin[i + 1][lt][ln] = roundz(Tmin[i][lt][ln] - 273.15,1);                Windspeed[i + 1][lt][ln] = roundz(Windspeed[i][lt][ln],1);                Radiation[i + 1][lt][ln]  = 1000 * roundz(86.400 * Radiation[i][lt][ln],1); // from W/m2 to J/m2/day                                /* Calculate the saturated vapour pressure */                Svap_Tmax = 6.108 * exp((17.27 * Tmax[i + 1][lt][ln]) / (237.3 + Tmax[i + 1][lt][ln]));                Svap_Tmin = 6.108 * exp((17.27 * Tmin[i + 1][lt][ln]) / (237.3 + Tmin[i + 1][lt][ln]));                Svap = (Svap_Tmax + Svap_Tmin) / 2.;                                /* Vapour pressure */                Vapour[i + 1][lt][ln] = roundz(0.01 * Vapour[i][lt][ln] * Svap,2); //hPa                                    //temp = 0.5 * (Tmax[i + 1][lt][ln] + Tmin[i + 1][lt][ln]);                //Vapour[i + 1][lt][ln]  = roundz(6.108 * exp((17.27*temp)/(237.3 + temp)),2); //hPa                                /* Set some limitations similar to PCSE */                if (Radiation[i + 1][lt][ln] > 39999999.) Radiation[i + 1][lt][ln] = 39999999.;                if (Tmax[i + 1][lt][ln] > 50.)    Tmax[i + 1][lt][ln] = 50.;                if (Rain[i + 1][lt][ln] > 249.5)  Rain[i + 1][lt][ln] = 249.5;            }        }    } return 1;}