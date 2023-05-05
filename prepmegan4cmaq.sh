#!/bin/bash
#------------------------------------------------
#author: Ramiro A. Espada. April, 2023.
#------------------------------------------------
export LC_NUMERIC="en_US.UTF-8"

#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

  srsInp="epsg:4326"	#(latlon spatial reference system of input files)

GRIDDESC_file="/home/usuario/runs/papila2019/cmaq/mcip/GRIDDESC" #"./GRIDDESC_example"
GRIDNAME="PAPILAGRID" #"Argentina"
   
#Input Files:
       laiv_files='./input/laiv2003'			 #LAIv=LAI/VegCover  (-) -  netCDFfiles path PREFIX

   crop_frac_file='./input/GF3aCrop.nc' 		 #frac of crop   cover (%) - netCDF file path
  grass_frac_file='./input/GF3aGrass.nc'		 #frac of grass  cover (%) - netCDF file path
  shrub_frac_file='./input/GF3aShrub.nc'		 #frac of shrubs cover (%) - netCDF file path
   tree_frac_file='./input/GF3aTree.nc'		         #frac of trees  cover (%) - netCDF file path

nl_tree_frac_file='./input/NTfrac_reorder_lat.nc'	 #frac of needleleaf trees (%) - netCDF file path
tp_tree_frac_file='./input/tropfrac_reorder_lat.nc'	 #frac of tropical   trees (%) - netCDF file path

     ecotype_file='./input/EVT3b.nc'			 #Gridded ecotypes ids         - netCDF file_path

     GtEcoEFfile="./db/GtEFbyEcotype.csv"		 #Emission factor of each GT grouped by Ecotype

#Input files for BDSNP:
arid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_arid.nc'	         #     arid mask (0/1)- netCDF file path
nonarid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_non_arid.nc'   # non arid mask (0/1)- netCDF file path
  landtype_files='./input/MEGAN31_Prep_Input_soil_191022/soil_landtype_'         # landtype           - netCDF files PREFIX
  nitro_depo_files='input/MEGAN31_Prep_Input_soil_191022/soil_nitrogen_mon'      # soil-NO deposition of each month (kg/m2/s) - netCDF files PREFIX
  #     fert_files='input/MEGAN31_Prep_Input_soil_191022/soil_fert_'             #Reservoir of N associated w/ manure and fertilizer (mg/m3) - netCDF files PREFIX

#-------------------------------------------------
#(0) Get grid & proj parameters from GRIDDESC:
read projName xorig yorig dx dy nx ny nz <<< $( sed -n "/${GRIDNAME}/{n;p;q}" "$GRIDDESC_file" )
read COORDTYPE P_ALP P_BET P_GAM XCENT YCENT <<< $( sed -n "/${projName}/{n;p;q}" "$GRIDDESC_file" )

  if [ $COORDTYPE == 1 ]; then     #Geographic:
   srsOut="+proj=latlong +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 2 ]; then     #Lambert Conformal Conic:
   srsOut="+proj=lcc +lat_1=$P_ALP +lat_2=$P_BET +lon_0=$P_GAM +lat_0=$YCENT +a=6370000.0 +b=6370000.0 +units=m"
elif [ $COORDTYPE == 3 ]; then     #General Mercator
   echo  "proyección: 3 (General Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 4 ]; then     #General tangent Stereografic
   srsOut="+proj=stere +lat_ts=$P_ALP +lat_0=$P_BET +lon_0=$P_GAM +a=6370000.0 +b=6370000.0 +k_0=1.0 +units=m"	#(!) Casi seguro que está mal
elif [ $COORDTYPE == 5 ]; then     #UTM
   echo  "proyección: 5 (Universal Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 6 ]; then     #Polar Secant Stereographic
   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0 +units=m"
elif [ $COORDTYPE == 7 ]; then     #Equatorial Mercator
   srsOut="+proj=merc +lat_ts=$P_ALP +lon_0=$P_GAM +a=6370000.0 +b=6370000.0 +units=m" 
elif [ $COORDTYPE == 8 ]; then     #Transverse Mercator
   echo  "proyección: 8 (Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 9 ]; then     #Lambert Azimuthal Equal-Area
   echo  "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada en esta aplicación."; stop
else
   echo  "codigo de proyección invalido. COORDTYPE: $COORDTYPE"; stop
fi;

echo "SRS of Input  Grids: $srsInp"
echo "SRS of Output Grids: $srsOut"
#Grilla
#xmin=$xorig;ymin=$yorig; xmax=$( bc -l <<<"${xorig}*-1 " );ymax=$( bc -l <<<"${yorig}*-1 " )
xmin=$xorig;ymin=$yorig; xmax=$( bc -l <<<"${xorig}+$nx*$dx " );ymax=$( bc -l <<<"${yorig}+$ny*$dy " )
#------------------------------------------------
#(1) Re-gridding: Agarrar los inputs y regrillarlos (e interpolar) segun GRIDDESC.
#interpolation methods: near (default), bilinear, cubic, cubicspline, lanczos, average, rms, mode,  max, min, med, Q1, Q3, sum
if [ ! -d "./tmp_grids" ]; then
    mkdir tmp_grids
fi

echo "Regridding input files..."

echo " $crop_frac_file    -> crop.nc   "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $crop_frac_file ./tmp_grids/crop.nc 
echo " $grass_frac_file   -> grass.nc  "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $grass_frac_file ./tmp_grids/grass.nc 
echo " $shrub_frac_file   -> shrub.nc  "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $shrub_frac_file ./tmp_grids/shrub.nc
echo " $tree_frac_file    -> tree.nc   "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $tree_frac_file ./tmp_grids/tree.nc 
echo " $nl_tree_frac_file -> nl_tree.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $nl_tree_frac_file ./tmp_grids/nl_tree.nc 
echo " $tp_tree_frac_file -> tp_tree.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $tp_tree_frac_file ./tmp_grids/tp_tree.nc
echo " $ecotype_file      -> ecotype.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode     -f "NetCDF"  $ecotype_file ./tmp_grids/ecotype.nc 

echo " $laiv_files        -> laiv.nc   "
for MM in $(seq --format='%02.0f' 1 1 12)
do
	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" ${laiv_files}${MM}* ./tmp_grids/laiv${MM}.nc 
done

echo "Regridding input files for BDSNP..."

echo " $arid_file     ->    arid.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF"  $arid_file     ./tmp_grids/arid.nc 
echo " $nonarid_file  -> nonarid.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF"  $nonarid_file  ./tmp_grids/nonarid.nc 

#Para Landtype hay 24 archivos(creo que son 24clases de land, y hay que ir sumandolos)
echo "   $landtype_files        -> landtype.nc   "
for LT in $(seq --format='%02.0f' 1 1 24)
do
	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF" ${landtype_files}${LT}* ./tmp_grids/landtype${LT}.nc 
done

echo " $nitro_depo_files        -> nitro.nc   "
for MM in $(seq --format='%02.0f' 1 1 12)
do
	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" ${nitro_depo_files}${MM}* ./tmp_grids/nitro${MM}.nc 
done

#(2) Llamar a programa en Fortran que ensamble todas las grillas y hace los calculos pertinentes
#Armo namelist input para prepmegan4cmaq.exe
cat << EOF > example.inp 
&control
start_date='${start_date}',
  end_date='${end_date}',

    griddesc_file='${GRIDDESC_file}'
         gridname='${GRIDNAME}'
        
        laiv_files='./tmp_grids/laiv',				

    crop_frac_file='./tmp_grids/crop.nc', 
   grass_frac_file='./tmp_grids/grass.nc',
   shrub_frac_file='./tmp_grids/shrub.nc',
    tree_frac_file='./tmp_grids/tree.nc', 
 nl_tree_frac_file='./tmp_grids/nl_tree.nc',		
 tp_tree_frac_file='./tmp_grids/tp_tree.nc',		
 bl_tree_frac_file='./tmp_grids/bl_tree.nc',		

      ecotype_file='./tmp_grids/ecotype.nc',				

      GtEcoEF_file='${GtEcoEFfile}',
!BDSNP
         arid_file='./tmp_grids/arid.nc',     
        narid_file='./tmp_grids/nonarid.nc',  
           lt_file='./tmp_grids/landtype',  

       nitro_files='./tmp_grids/nitro',
       !fert_files='./tmp_grids/fert'
/
EOF

./prepmegan4cmaq.exe < example.inp


#####################################################################
#GRIDDESC proj params (para usarlo para armar los srsOut):
#
#GT==2:
#	Lambert Conformal Conic 
#	LAMGRD3 = 2
#	PROJ_ALPHA ≤ PROJ_BETA are the two latitudes which determine the projection cone; PROJ_GAMMA is the central meridian.
#(X_CENT,Y_CENT) are the (lon ,lat) coordinates for the center (0,0) of the Cartesian coordinate system.
#	Coordinate units are meters. GCTP projection 4.
#	
#GT== 4:
#	Stereographic 	
#	STEGRD3 = 4 	general tangent stereographic
#	PROJ_ALPHA is the "true latitude", the latitude at which the stereographic plane is secant to the Earth (or 90, if the projection is a North Polar tangent stereographic. PROJ_BETA is the angle of rotation of the Y-axis relative to the Greenwich meridian, i.e., the longitude meridian which becomes the negative Y axis. PROJ_GAMMA is unused. (X_CENT,Y_CENT) are the (lon ,lat) coordinates for the center (0,0) of the Cartesian coordinate system.
#	Coordinate units are meters.
#	GCTP projection 10.
#GT==6
##	Polar
##	POLGRD3 = 6
##	polar secant stereographic
##	PROJ_ALPHA is 1.0 for North Polar, -1.0 for South Polar. PROJ_BETA is the secant latitude (latitude of true scale), PROJ_GAMMA is the central meridian.(X_CENT,Y_CENT) are the (lon ,lat) coordinates for the center (0,0) of the Cartesian coordinate system.
##	Coordinate units are meters.
##	GCTP projection 6. 
#GT==7:
#	Equatorial Mercator
#	EQMGRD3 = 7
#	PROJ_ALPHA is the latitude of true scale. PROJ_BETA is unused.PROJ_GAMMA is the longitude of the central meridian.(X_CENT,Y_CENT) are the (lon ,lat) coordinates for the center (0,0) of the Cartesian coordinate system.
#	Coordinate units are meters.
#	GCTP projection 5.
#GT==8:
#	Transverse Mercator
#	TRMGRD3 = 8
#	PROJ_ALPHA is the latitude of true scale. PROJ_BETA is the scale factor at the central meridian;PROJ_GAMMA is the longitude of the central meridian. (X_CENT,Y_CENT) are the (lon ,lat) coordinates for the center (0,0) of the Cartesian coordinate system. 
#	Coordinate units are meters.
#	GCTP projection 9.
