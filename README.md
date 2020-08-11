# wofost_copernicus (netcdf)

run as follows: ./wofost list.txt meteolist.txt

Example list.txt:
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_1.dat  Site/France.site 02-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_2.dat  Site/France.site 02-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_3.dat  Site/France.site 02-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_4.dat  Site/France.site 03-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_5.dat  Site/France.site 03-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_6.dat  Site/France.site 03-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_7.dat  Site/France.site 04-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_8.dat  Site/France.site 04-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_9.dat  Site/France.site 04-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_10.dat  Site/France.site 05-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_11.dat  Site/France.site 05-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-1.cab Soil/ec3.new  Management/Fr_mg_12.dat  Site/France.site 05-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_1.dat  Site/France.site 02-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_2.dat  Site/France.site 02-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_3.dat  Site/France.site 02-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_4.dat  Site/France.site 03-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_5.dat  Site/France.site 03-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_6.dat  Site/France.site 03-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_7.dat  Site/France.site 04-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_8.dat  Site/France.site 04-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_9.dat  Site/France.site 04-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_10.dat  Site/France.site 05-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_11.dat  Site/France.site 05-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-2.cab Soil/ec3.new  Management/Fr_mg_12.dat  Site/France.site 05-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-3.cab Soil/ec3.new  Management/Fr_mg_1.dat  Site/France.site 02-05 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-3.cab Soil/ec3.new  Management/Fr_mg_2.dat  Site/France.site 02-15 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-3.cab Soil/ec3.new  Management/Fr_mg_3.dat  Site/France.site 02-25 0  ../Output/sb-G1.csv
../Parameters/ Crop/sb-G1-3.cab Soil/ec3.new  Management/Fr_mg_4.dat  Site/France.site 03-05 0  ../Output/sb-G1.csv

Eaxmple meteolist.txt:
../Forcing/France/ 1987 2016 ../Parameters/Mask/springbarley_area_rh_p.org
Tmin_daily_WFDEI_1987-2016.nc TMIN Tmin
Tmax_daily_WFDEI_1987-2016.nc TMAX  Tmax
SWdown_daily_WFDEI_1987-2016.nc RADIATION SWdown
Rainf_daily_WFDEI_1987-2016.nc RAIN Rainf
Wind_daily_WFDEI_1987-2016.nc WINDSPEED Wind
Vap_daily_WFDEI_1987-2016.nc VAPOUR Vap
