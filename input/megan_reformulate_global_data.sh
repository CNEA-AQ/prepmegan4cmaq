#!/bin/bash

#
#Script para re-estructurar input files globales para preprocesador de MEGAN:
#

# General specs:
# - netCDF files: domain global, grilla regular, ~1km resol. 	
#	- lon = 43200 ; (-180 a 180 degrees )	orden creciente!
#	- lat = 16800 ;	(-57  a  82 degrees )   orden creciente!
#
#--------------------------------------------------------------------------------------------------------------------------------------------
# Original:			         Deseado:							
# ├── EVT3b.nc                           ├── veg_Ecotypes.nc	         [int  ] (ecotype)					       ( 2Gb)
# ├── GF3aCrop.nc                        ├── veg_GrowthFormFractions.nc  [byte ] (crop, grass, shrub, tree, nl_tree, bl_tree, tp_tree) ( 5Gb)
# ├── GF3aGrass.nc                       ├── veg_LAIv.nc                 [short] (laiv01,...,laiv12)				       (    )
# ├── GF3aShrub.nc                       ├── soil_Landtype.nc            [byte ] (landtype,arid, nonarid)			       (    )
# ├── GF3aTree.nc                        ├── soil_Nitro.nc               [float] (nitro01,...,nitro12)				       (    )
# ├── NTfrac_reorder_lat.nc              └── soil_Fert.nc                [float] (fert01,...,fert365)				       (    )
# ├── tropfrac_reorder_lat.nc            											--------------
# ├── laiv200301_30sec.nc                                                                                                       Total. (  Gb)
# ├── ...                                
# ├── laiv200312_30sec.nc                
# └── MEGAN31_Prep_Input_soil_191022     
#     ├── soil_climate_arid.nc           
#     ├── soil_climate_non_arid.nc       
#     ├── soil_fert_001.nc               
#     ├── soil_fert_002.nc               
#     ├── ...                            
#     ├── soil_fert_366.nc               
#     ├── soil_landtype_01.nc            
#     ├── ...                            
#     ├── soil_landtype_24.nc            
#     ├── soil_nitrogen_mon01.nc         
#     ├── ...                            
#     └── soil_nitrogen_mon12.nc         
#--------------------------------------------------------------------------------------------------------------------------------------------

#----------------------------
#Ecotypes. Está OK. Solo le cambio el nombre.
mv EVT3b.nc veg_Ecotypes.nc
#----------------------------
#Growtht Form:

#NTfrac esta invertido en la direccion-y (respecto de los demas) y es tipo "float", asi que lo convierto para que me quede como el resto.
ncpdq -a -lat NTfrac_reorder_lat.nc NTfrac.nc
ncap2 -4 -O -s "where(NTfrac < 0) NTfrac=0.0; where(NTfrac>0) NTfrac=NTfrac*100;" NTfrac.nc tmp.nc
ncap2 -4 -O -s "NTfrac = byte(NTfrac)" tmp.nc NTfrac.nc
#tropFrac es solo (0/1).

ncks -O -v "lat"       GF3aCrop.nc  tmp.nc
ncks -A -v "lon"       GF3aCrop.nc  tmp.nc
ncks -A -v "m20crop"   GF3aCrop.nc  tmp.nc
ncks -A -v "m202grass" GF3aGrass.nc tmp.nc
ncks -A -v "m204shrub" GF3aShrub.nc tmp.nc
ncks -A -v "m20tree"   GF3aTree.nc  tmp.nc
ncks -A -v "NTfrac"    NTfrac.nc    tmp.nc
ncks -A -v "tropfrac"  tropfrac.nc  tmp.nc

ncpdq -a -lat tmp.nc GrowthFormFracions.nc	#ordenar lat dimension (en orden creciente)

rm tmp.nc #limpio.
#----------------------------
# LAIv
# Hay que agruparlos todos, cambiar el missing value y luego convertir a [short]
ncks -4 -O -v "lat" laiv200301_30sec.nc tmp.nc
ncks -4 -A -v "lon" laiv200301_30sec.nc tmp.nc

NAMES=(" " "LAI_for_Jan_2003_\(m2_per_m2\)" "LAI_for_Feb_2003_\(m2_per_m2\)" "LAI_for_Mar_2003_\(m2_per_m2\)" "LAI_for_Apr_2003_\(m2_per_m2\)" "LAI_for_May_2003_\(m2_per_m2\)" "LAI_for_Jun_2003_\(m2_per_m2\)" "LAI_for_Jul_2003_\(m2_per_m2\)" "LAI_for_Aug_2003_\(m2_per_m2\)" "LAI_for_Sep_2003_\(m2_per_m2\)" "LAI_for_Oct_2003_\(m2_per_m2\)" "LAI_for_Nov_2003_\(m2_per_m2\)" "LAI_for_Dec_2003_\(m2_per_m2\)")
for mm in $(seq --format="%02.0f" 1 12)
do
	m="$(printf "%d" ${mm})";

	file=laiv2003${mm}_30sec.nc
	echo "$file, ${NAMES[$m]}"
	ncks -4 -A -v "LAI_*" $file tmp.nc
	ncrename -v ${NAMES[$m]},laiv${mm} tmp.nc

done

#----------------------------
# Landtype
# Hay que sumarlos a todos, y luego convertir a [byte].
# Luego agregar las dos capas "arid" y "nonarid"
ncks -4 -O -v "lat" MEGAN31_Prep_Input_soil_191022/soil_landtype_01.nc tmp.nc
ncks -4 -A -v "lon" MEGAN31_Prep_Input_soil_191022/soil_landtype_01.nc tmp.nc

for lt0 in $(seq --format="%02.0f" 1 24)
do
	lt="$(printf "%d" ${lt0})";

	file=MEGAN31_Prep_Input_soil_191022/soil_landtype_${lt0}.nc
        ncks -A -v LANDFRAC $file -o tmp.nc
        #ncap2 -h -A -C -s "landtype=landtype+LANDFRAC*${lt};LANDFRAC=0.0" tmp.nc
	#ncap2 -O -s 'landtype01=file1_landtype[:]*1;' -s 'landtype02=file2_landtype[:]*2;' ... -s 'finalLandtype=landtype01+landtype02+...;' output.nc
done
			
                        



#----------------------------
# Nitro



#----------------------------
# Ferti




