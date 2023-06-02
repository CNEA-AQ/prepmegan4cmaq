#!/bin/bash

#
#Script para re-estructurar input files globales para preprocesador de MEGAN:
#

# General specs:
# - netCDF files: domain global, grilla regular, ~1km resol. 	
#	- lon = 43200 ; (-180 a 180 degrees )	orden creciente!
#	- lat = 16800 ;	(-57  a  82 degrees )   orden creciente!
#
#Data types: 
#
# byte (          -128 to 127            )
# short(       -32,768 to 32,767         )
# int  (-2,147,483,648 to 2,147,483,647  )
# float(                                 )
#--------------------------------------------------------------------------------------------------------------------------------------------
# Original:			                  Deseado:                                                                                  
# ├── EVT3b.nc              [btye  ]  ( 1.5Gb)    ├── veg_Ecotypes.nc   [int  ] (ecotype)				              (1.5Gb)
# ├── GF3aCrop.nc           [btye  ]  ( 0.7Gb)    ├── veg_GrowthForm.nc [byte ] (crop,grass,shrub,tree,nl_tree,bl_tree,tp_tree)       (  3Gb)
# ├── GF3aGrass.nc          [btye  ]  ( 0.7Gb)    ├── veg_LAIv.nc  *    [short] (laiv01,...,laiv12)		                      (17.Gb)
# ├── GF3aShrub.nc          [btye  ]  ( 0.7Gb)    ├── soil_Climate.nc   [byte ] (arid, nonarid)                                       (     )
# ├── GF3aTree.nc           [btye  ]  ( 0.7Gb)    ├── soil_Landtype.nc  [byte ] (landtype)		                              (  4Mb)
# ├── NTfrac_reorder_lat.nc [float ]  ( 2.0Gb)    ├── soil_Nitro.nc*    [float] (nitro01,...,nitro12)		                      (     )
# ├── tropfrac_reorder_lat.nc [byte]  ( 0.6Gb)    └── soil_Ferti.nc*    [float] (fert01,...,fert365)		                      (     )
# ├── laiv200301_30sec.nc      
# ├── ...                   [float ]  (      )
# ├── laiv200312_30sec.nc                    
# └── MEGAN31_Prep_Input_soil_191022         
#     ├── soil_climate_arid.nc        (  1Mb )
#     ├── soil_climate_non_arid.nc    (  1Mb )
#     ├── soil_fert_001.nc                   
#     ├── soil_fert_002.nc                   
#     ├── ...                         ( 366Mb)
#     ├── soil_fert_366.nc                   
#     ├── soil_landtype_01.nc        
#     ├── ...                         (  24Mb)         
#     ├── soil_landtype_24.nc         
#     ├── soil_nitrogen_mon01.nc         
#     ├── ...                         (  12Mb)
#     └── soil_nitrogen_mon12.nc         
#                                    ----------                                                                                --------------
#                              Total. ( ~8 Gb)                                                                                  Total.( ~6Gb)
#--------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------
#Ecotypes. Está OK. Solo le cambio el nombre.
#mv EVT3b.nc veg_Ecotypes.nc
ncpdq -4 -O -a -lat EVT3b.nc veg_Ecotypes.nc
#ncap2   -h -O -s "where(ecotype == 32767) ecotype=-999; where(ecotype < 0) ecotype=-999;" veg_Ecotypes.nc veg_Ecotypes.nc
#ncatted -h -O -a missing_value,ecotype,o,s,-999 veg_Ecotypes.nc
ncatted -h -O -a esri_pe_string,ecotype,d,s,"" veg_Ecotypes.nc
ncatted -h -O -a units,ecotype,o,c,"" veg_Ecotypes.nc
ncatted -h -O -a Conventions,global,d,c,"" veg_Ecotypes.nc
ncatted -h -O -a Source_Software,global,d,c,"" veg_Ecotypes.nc

#----------------------------
#Growtht Form:
#NTfrac esta invertido en la direccion-y (respecto de los demas) y es tipo "float", asi que lo convierto para que me quede como el resto.
ncpdq -4 -O -a -lat NTfrac_reorder_lat.nc NTfrac.nc
ncap2 -4 -O -s "where(NTfrac < 0) NTfrac=-128; where(NTfrac > 999) NTfrac=-128; where(NTfrac>0) NTfrac=NTfrac*100; NTfrac = byte(NTfrac)" NTfrac.nc tmp.nc
ncap2 -4 -O -s "NTfrac = byte(NTfrac)" tmp.nc NTfrac.nc

ncks -h -O -v "lat"       GF3aCrop.nc  tmp.nc
ncks -h -A -v "lon"       GF3aCrop.nc  tmp.nc
ncks -h -A -v "m20crop"   GF3aCrop.nc  tmp.nc
ncks -h -A -v "m202grass" GF3aGrass.nc tmp.nc
ncks -h -A -v "m204shrub" GF3aShrub.nc tmp.nc
ncks -h -A -v "m20tree"   GF3aTree.nc  tmp.nc
ncks -h -A -v "NTfrac"    NTfrac.nc    tmp.nc
#ncks -h -A -v "tropfrac"  tropfrac.nc  tmp.nc
ncks -h -A -v "tropfrac"  tropfrac_reorder_lat.nc  tmp.nc

ncap2 -s 'tropfrac = 100*tropfrac' -O tmp.nc tmp.nc		#tropfrac [0/1] -> (0,100)
#Ajuste de metadata:
ncatted -h -O -a long_name,m20crop,d,c,"crop cover fraction" tmp.nc
ncatted -h -O -a long_name,m202grass,d,c,"grass cover fraction" tmp.nc
ncatted -h -O -a long_name,m204shrub,d,c,"shrub cover fraction" tmp.nc
ncatted -h -O -a long_name,m20tree,d,c,"tree cover fraction" tmp.nc
ncatted -h -O -a long_name,NTfrac,d,c,"needle leaf tree cover fraction" tmp.nc
ncatted -h -O -a long_name,tropfrac,d,c,"tropical tree cover fraction" tmp.nc

ncatted -h -O -a esri_pe_string,m20crop,d,s,"" tmp.nc
ncatted -h -O -a esri_pe_string,m202grass,d,s,"" tmp.nc
ncatted -h -O -a esri_pe_string,m204shrub,d,s,"" tmp.nc
ncatted -h -O -a esri_pe_string,m20tree,d,s,"" tmp.nc
ncatted -h -O -a esri_pe_string,NTfrac,d,s,"" tmp.nc
ncatted -h -O -a esri_pe_string,tropfrac,d,s,"" tmp.nc

ncatted -h -O -a missing_value,tropfrac,o,b,-128 tmp.nc
ncatted -h -O -a missing_value,NTfrac,o,b,-128 tmp.nc
ncatted -h -O -a units,NTfrac,o,c,"percent" tmp.nc

ncatted -h -O -a Source_Software,global,d,c,"" tmp.nc
ncatted -h -O -a Conventions,global,d,c,"" tmp.nc
ncatted -h -O -a NCO,global,d,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc

ncpdq -h -a -lat tmp.nc GrowthFormFracions.nc	#ordenar lat dimension (en orden creciente)

ncrename -h -O -v m20crop,crop    GrowthFormFracions.nc
ncrename -h -O -v m202grass,grass GrowthFormFracions.nc
ncrename -h -O -v m204shrub,shrub GrowthFormFracions.nc
ncrename -h -O -v m20tree,tree    GrowthFormFracions.nc
ncrename -h -O -v NTfrac,nl_tree  GrowthFormFracions.nc
ncrename -h -O -v tropfrac,trop_tree GrowthFormFracions.nc

rm tmp.nc #limpio.
#----------------------------
# LAIv
# Hay que agruparlos todos, cambiar el missing value y luego convertir a [short]
ncks -4 -O -v "lat" laiv200301_30sec.nc tmp.nc
ncks -4 -A -v "lon" laiv200301_30sec.nc tmp.nc

NAMES=(" " "LAI_for_Jan_2003_(m2_per_m2)" "LAI_for_Feb_2003_(m2_per_m2)" "LAI_for_Mar_2003_(m2_per_m2)" "LAI_for_Apr_2003_(m2_per_m2)" "LAI_for_May_2003_(m2_per_m2)" "LAI_for_Jun_2003_(m2_per_m2)" "LAI_for_Jul_2003_(m2_per_m2)" "LAI_for_Aug_2003_(m2_per_m2)" "LAI_for_Sep_2003_(m2_per_m2)" "LAI_for_Oct_2003_(m2_per_m2)" "LAI_for_Nov_2003_(m2_per_m2)" "LAI_for_Dec_2003_(m2_per_m2)")
for mm in $(seq --format="%02.0f" 1 12)
do
	m=${mm#0} #"$(printf "%d" ${mm})";
	file=laiv2003${mm}_30sec.nc
	echo "$file, ${NAMES[$m]}"
	ncks -4 -A -v "LAI_*" $file tmp.nc
	ncrename -v ${NAMES[$m]},laiv${mm} tmp.nc
        ncatted -h -O -a esri_pe_string,laiv${mm},d,c,"" tmp.nc
        ncatted -h -O -a units,laiv${mm},o,c,"m2_per_m2" tmp.nc
        #ncatted -h -O -a missing_value,laiv${mm},o,s,"-32768" tmp.nc
done
ncatted -h -O -a Source_Software,global,d,c,"" tmp.nc
ncatted -h -O -a Conventions,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc

mv tmp.nc veg_LAIv.nc
#----------------------------
# Landtype
# Hay que sumarlos a todos, y luego convertir a [byte].
ncks -O -v "lat" -x -C MEGAN31_Prep_Input_soil_191022/soil_landtype_01.nc tmp.nc
ncks -A -v "lon" -x -C MEGAN31_Prep_Input_soil_191022/soil_landtype_01.nc tmp.nc

ncks -A -v LANDFRAC -x -C MEGAN31_Prep_Input_soil_191022/soil_landtype_01.nc tmp.nc
ncap2 -A -s "landtype=(LANDFRAC)" tmp.nc tmp.nc

for lt0 in $(seq --format="%02.0f" 2 24)
do
	file=MEGAN31_Prep_Input_soil_191022/soil_landtype_${lt0}.nc
	lt=${lt0#0}
        ncks  -A -v LANDFRAC $file tmp.nc
	ncap2 -A -s "landtype+=LANDFRAC*${lt}" tmp.nc tmp.nc
done
ncap2 -h -A -C -s "landtype=byte(landtype)" tmp.nc tmp.nc
ncatted -h -O -a Source_Software,global,d,c,"" tmp.nc
ncatted -h -O -a Conventions,global,d,c,"" tmp.nc
ncatted -h -O -a NCO,global,d,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc
ncatted -h -O -a title,global,o,c,"Landtypes for MODIS-KOPPEN Biome Types" tmp.nc

ncatted -h -O -a long_name,landtype,o,c,"Dominant Landtype categories for MODIS-KOPPEN Biome Types" tmp.nc

ncks -O -x -v LANDFRAC tmp.nc soil_landtype.nc
rm tmp.nc
#----------------------------
# Climate
# Luego agregar las dos capas "arid" y "nonarid"

ncks  -O -v arid    MEGAN31_Prep_Input_soil_191022/soil_climate_arid.nc     tmp.nc
ncks  -A -v non_arid MEGAN31_Prep_Input_soil_191022/soil_climate_non_arid.nc tmp.nc

ncap2 -4 -O -s "arid = byte(arid)" tmp.nc tmp.nc
ncap2 -4 -O -s "non_arid = byte(non_arid)" tmp.nc tmp.nc

ncatted -h -O -a conventions,global,d,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc

mv tmp.nc soil_climate.nc

#----------------------------
# Nitro
rm tmp.nc
for mm in $(seq --format="%02.0f" 1 12)
do
	m=${mm#0}
	month=$(date -d "2023-${mm}-01" +'%b')

	file=MEGAN31_Prep_Input_soil_191022/soil_nitrogen_mon${mm}.nc 
	echo "$file, nitrogen"

	ncks -A -v "nitrogen" $file tmp.nc
	ncrename -v "nitrogen",nitro${mm} tmp.nc
	
	ncatted -h -O -a long_name,nitro${mm},o,c,"Monthly Dry and wet deposition of N (${month})" tmp.nc
done
ncatted -h -O -a conventions,global,d,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc

mv tmp.nc soil_nitro.nc
#----------------------------
# Ferti
for ddd in $(seq --format="%03.0f" 1 366)
do
	file=MEGAN31_Prep_Input_soil_191022/soil_fert_${ddd}.nc 
	echo "$file, fert"
	ncks -h -4 -A -v "fert" $file tmp.nc
	ncrename -h -v "fert",fert${ddd} tmp.nc

done
ncatted -h -O -a conventions,global,d,c,"" tmp.nc
ncatted -h -O -a history_of_appended_files,global,d,c,"" tmp.nc
ncatted -h -O -a history,global,o,c,"" tmp.nc

mv tmp.nc soil_fert.nc
#============================
#Probblems i have faced:
# + 
# + 
# + 
# + 
 

