#!/bin/bash

case $1 in
(-download)
	#dates=("2009-03-16" "2009-03-28" "2009-04-05" )
        dates=("2009-08-05" "2010-07-11" "2010-07-12" "2015-04-23" "2016-11-03" "2017-09-13")
	date_ini="2009-03-28"	#YYYY-MM-DD
	date_end="2009-03-16"	#YYYY-MM-DD
	cc=61			#collection number
	myToken="ZXNwYWRhOlpYTndZV1JoUUdGbmNtOHVkV0poTG1GeToxNjIxMjgyMzA4OjgxYWY4Y2Q0ZTBjODNiOGEzMjdjODBiMDRlMjlmOTFmZjE3MmE3Yjg"
	for date in ${dates[@]}
	do
		read yr mo day jdy<<<`date -u -d"${date} 00 hours" +"%Y %m %d %j"`
		
	wget -e robots=off -m -np -nH -A "MYD04_L2.A$yr$jdy.1[678]*.hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/$cc/MYD04_L2/$yr/$jdy/" --header "Authorization: Bearer ${myToken}" -P .
	#wget -e robots=off -m -np -nH -A "MOD04_L2.A$yr$jdy.19*.hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/$cc/MOD04_L2/$yr/$jdy/" --header "Authorization: Bearer ${myToken}" -P .
	done;
;;
(-extract)

#paths:
paths=($(ls -d archive/allData/61/MYD04_L2/20*/*))

for path in ${paths[@]}
do
	files=($(ls $path))
	for file in ${files[@]}
	do
		product=Mass_Concentration_Land
		#product=Mass_Concentration_Ocean
		#product=Aerosol_Type_Land
		#product=Optical_Depth_Land_And_Ocean
		#product=Image_Optical_Depth_Land_And_Ocean
		#product=Deep_Blue_Aerosol_Optical_Depth_550_Land
		#product=Deep_Blue_Angstrom_Exponent_Land
		#product=Deep_Blue_Spectral_Aerosol_Optical_Depth_Land
		#product=Deep_Blue_Cloud_Fraction_Land
		#product=Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate
		#product=AOD_550_Dark_Target_Deep_Blue_Combined

		gdal_translate -ot FLOAT64 -of GTiff HDF4_EOS:EOS_SWATH:"${path}/${file}":mod04:${product} tmp/${file}.tif
	done
	
	#Unir todos los tiff y cropearlos (falta esto ultimo)
	xmin=-90; xmax=-50
       	ymin=-80; ymax=-20
	gdalwarp -ot FLOAT64 -te $xmin $ymin $xmax $ymax tmp/* tmp.tiff
	#Corregir (escalar) valores
	name=$(gdalinfo tmp.tiff | grep "RANGEBEGINNINGDATE=" | sed 's/^.*=//')
        scale=$(gdalinfo tmp.tiff | grep "scale_factor=" | sed 's/^.*=//')
	gdal_calc.py -A tmp.tiff --type='Float64' --outfile=${product}_${name}.tiff --calc="A*${scale}"
	
	rm tmp.tiff tmp/*
done
;;
esac
# gdalinfo HDF4_EOS:EOS_SWATH:"archive/allData/61/MYD04_L2/2009/075/MYD04_L2.A2009075.1805.061.2018042105351.hdf":mod04:Image_Optical_Depth_Land_And_Ocean | grep BOUND
#gdal_merge.py -n 9999.990234375 -of GTiff -o /path/to/myoutputfile.tif /path/to/inputfile1.tif /path/to/inputfile2.tif /path/to/inputfileX.tif
#gdalwarp -s_srs EPSG:4326 -t_srs epsg:4326 -dstnodata 0 -of GTiff /path/to/myoutputfile.tif /path/to/newoutputfile.tif


