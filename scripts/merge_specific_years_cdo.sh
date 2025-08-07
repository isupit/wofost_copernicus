#!/bin/bash
# A script to merge specific yearly NetCDF files (1986-2015) for multiple variables
# using CDO (Climate Data Operators). Uses multi-step approach to avoid chaining issues.
#
# This version:
# 1. Merges files.
# 2. Applies unit conversions, preserving time variable.
# 3. Adjusts time to seconds since 1970-01-01.
# 4. Handles fill values to prevent HDF5 errors.
# 5. Sets metadata to match China dataset, step-by-step.
echo "--- Starting CDO merge and process workflow for 1986-2015 ---"
# Define common variables
TARGET_TIME_UNITS="seconds since 1970-01-01 00:00:00"
TARGET_CALENDAR="standard"
SECONDS_OFFSET=$((5844 * 86400)) # Seconds from 1970-01-01 to 1986-01-01
# --- 1. SWdown (Shortwave Downward Radiation) ---
echo ""
echo "Processing: SWdown"
# Step 1a: MERGE
cdo -O mergetime SWdown_daily_WFDEI_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown merge."; exit 1; fi
# Step 1b: CHECK time variable after merge
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 1c: CONVERT units (W/m2 to KJ m-2 day-1)
cdo -O mulc,86.4 temp_merged.nc temp_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown units."; exit 1; fi
# Step 1d: CHECK time variable after units
echo "Checking time variable in temp_units.nc"
ncdump -h temp_units.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_units.nc"; exit 1; }
# Step 1e: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_units.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown time."; exit 1; fi
# Step 1f: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown setrtomiss."; exit 1; fi
cdo -O setmissval,0 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown setmissval."; exit 1; fi
# Step 1g: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown calendar."; exit 1; fi
# Step 1h: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown time units."; exit 1; fi
# Step 1i: SET SWdown metadata
cdo -O -f nc -setattribute,SWdown@long_name="Rate of surface downward shortwave radiation (rate over the previous 24 hours)" -setunit,"KJ m-2 day-1" temp_time_units.nc SWdown_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on SWdown metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_units.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc
echo "SUCCESS: Created SWdown_daily_1986-2015.nc"
# --- 2. Rainf (Rainfall) ---
echo ""
echo "Processing: Rainf"
# Step 2a: MERGE
cdo -O mergetime Rainf_daily_WFDEI_CRU_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf merge."; exit 1; fi
# Step 2b: CHECK time variable
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 2c: CONVERT units (kg/m2/s to mm)
cdo -O mulc,86400 temp_merged.nc temp_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf units."; exit 1; fi
# Step 2d: CHECK time variable
echo "Checking time variable in temp_units.nc"
ncdump -h temp_units.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_units.nc"; exit 1; }
# Step 2e: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_units.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf time."; exit 1; fi
# Step 2f: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf setrtomiss."; exit 1; fi
cdo -O setmissval,0 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf setmissval."; exit 1; fi
# Step 2g: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf calendar."; exit 1; fi
# Step 2h: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf time units."; exit 1; fi
# Step 2i: SET Rainf metadata
cdo -O -f nc -setattribute,Rainf@long_name="Sum of surface precipitation (sum over the previous 24 hours)" -setunit,"mm" temp_time_units.nc Rainf_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Rainf metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_units.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc
echo "SUCCESS: Created Rainf_daily_1986-2015.nc"
# --- 3. Tmax (Maximum Temperature) ---
echo ""
echo "Processing: Tmax"
# Step 3a: MERGE
cdo -O mergetime Tmax_daily_WFDEI_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax merge."; exit 1; fi
# Step 3b: CHECK time variable
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 3c: CONVERT units (K to degree C)
cdo -O -expr,'Tmax=Tmax-273.15' temp_merged.nc temp_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax units."; exit 1; fi
# Step 3d: CHECK time variable
echo "Checking time variable in temp_units.nc"
ncdump -h temp_units.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_units.nc"; exit 1; }
# Step 3e: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_units.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax time."; exit 1; fi
# Step 3f: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax setrtomiss."; exit 1; fi
cdo -O setmissval,-9999 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax setmissval."; exit 1; fi
# Step 3g: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax calendar."; exit 1; fi
# Step 3h: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax time units."; exit 1; fi
# Step 3i: REMOVE bounds variable if present (using selname to exclude time_bnds)
cdo -O selname,Tmax,time,lat,lon temp_time_units.nc temp_bounds.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax bounds."; exit 1; fi
# Step 3j: SET Tmax metadata
cdo -O -f nc -setattribute,Tmax@long_name="Daily Maximum Near-Surface Air Temperature" -setunit,"degree C" temp_bounds.nc Tmax_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmax metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_units.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc temp_bounds.nc
echo "SUCCESS: Created Tmax_daily_1986-2015.nc"
# --- 4. Vap (Vapour Pressure) ---
echo ""
echo "Processing: Vap (from Vapour) using the robust method..."
# Step 4a: MERGE
cdo -O mergetime Vapour_daily_WFDEI_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap merge."; exit 1; fi
# Step 4b: CHECK time variable
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 4c: CALCULATE Pa->kPa and rename
cdo -O -expr,'Vap=Vapour/1000' temp_merged.nc temp_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap units."; exit 1; fi
# Step 4d: CHECK time variable
echo "Checking time variable in temp_units.nc"
ncdump -h temp_units.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_units.nc"; exit 1; }
# Step 4e: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_units.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap time."; exit 1; fi
# Step 4f: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap setrtomiss."; exit 1; fi
cdo -O setmissval,0 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap setmissval."; exit 1; fi
# Step 4g: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap calendar."; exit 1; fi
# Step 4h: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap time units."; exit 1; fi
# Step 4i: SET Vap metadata
cdo -O -f nc -setattribute,Vap@long_name="Instantaneous near surface vapour pressure (average over the previous 24 hours)" -setunit,"kPa" temp_time_units.nc Vap_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Vap metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_units.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc
echo "SUCCESS: Created Vap_daily_1986-2015.nc"
# --- 5. Tmin (Minimum Temperature) ---
echo ""
echo "Processing: Tmin"
# Step 5a: MERGE
cdo -O mergetime Tmin_daily_WFDEI_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin merge."; exit 1; fi
# Step 5b: CHECK time variable
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 5c: CONVERT units (K to degree C)
cdo -O -expr,'Tmin=Tmin-273.15' temp_merged.nc temp_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin units."; exit 1; fi
# Step 5d: CHECK time variable
echo "Checking time variable in temp_units.nc"
ncdump -h temp_units.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_units.nc"; exit 1; }
# Step 5e: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_units.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin time."; exit 1; fi
# Step 5f: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin setrtomiss."; exit 1; fi
cdo -O setmissval,-9999 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin setmissval."; exit 1; fi
# Step 5g: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin calendar."; exit 1; fi
# Step 5h: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin time units."; exit 1; fi
# Step 5i: REMOVE bounds variable if present (using selname to exclude time_bnds)
cdo -O selname,Tmin,time,lat,lon temp_time_units.nc temp_bounds.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin bounds."; exit 1; fi
# Step 5j: SET Tmin metadata
cdo -O -f nc -setattribute,Tmin@long_name="Daily Minimum Near-Surface Air Temperature" -setunit,"degree C" temp_bounds.nc Tmin_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Tmin metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_units.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc temp_bounds.nc
echo "SUCCESS: Created Tmin_daily_1986-2015.nc"
# --- 6. Wind (Wind Speed) ---
echo ""
echo "Processing: Wind"
# Step 6a: MERGE
cdo -O mergetime Wind_daily_WFDEI_{1986..2015}.nc temp_merged.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind merge."; exit 1; fi
# Step 6b: CHECK time variable
echo "Checking time variable in temp_merged.nc"
ncdump -h temp_merged.nc | grep "time" || { echo "FATAL ERROR: time variable not found in temp_merged.nc"; exit 1; }
# Step 6c: ADJUST time
cdo -O shifttime,${SECONDS_OFFSET}sec temp_merged.nc temp_time.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind time."; exit 1; fi
# Step 6d: REPLACE fill values
cdo -O setrtomiss,1e20,1e20 temp_time.nc temp_clean1.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind setrtomiss."; exit 1; fi
cdo -O setmissval,0 temp_clean1.nc temp_clean.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind setmissval."; exit 1; fi
# Step 6e: SET calendar
cdo -O -setcalendar,standard temp_clean.nc temp_cal.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind calendar."; exit 1; fi
# Step 6f: SET time units
cdo -O -setattribute,time@units="seconds since 1970-01-01 00:00:00" temp_cal.nc temp_time_units.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind time units."; exit 1; fi
# Step 6g: SET Wind metadata
cdo -O -f nc -setattribute,Wind@long_name="Instantaneous near surface wind speed (average over the previous 24 hours)" -setunit,"m s-1" temp_time_units.nc Wind_daily_1986-2015.nc
if [ $? -ne 0 ]; then echo "FATAL ERROR on Wind metadata."; exit 1; fi
# Cleanup
rm temp_merged.nc temp_time.nc temp_clean1.nc temp_clean.nc temp_cal.nc temp_time_units.nc
echo "SUCCESS: Created Wind_daily_1986-2015.nc"
echo ""
echo "--- All processing complete. ---"
