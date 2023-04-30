#!/bin/bash
#------------------------------------------------
#author: Ramiro A. Espada. April, 2023.
#------------------------------------------------
export LC_NUMERIC="en_US.UTF-8"

#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

   srsInp="epsg:4326"	#Cambiar <solo si> los archivos NO vienen en latlon.

GRIDDESC_file="./GRIDDESC_example"
GRIDNAME="Argentina"

#Input Files:
      laiv_files='./input/laiv2003'			#(LAIv=LAI/VEGCOVER) file path just the name before the month indicator!

      crop_frac_file='./input/GF3aCrop.nc' 		#frac of crop cover       - file path
     grass_frac_file='./input/GF3aGrass.nc'		#frac of grass cover      - file path
     shrub_frac_file='./input/GF3aShrub.nc'		#frac of shrubs cover     - file path
      tree_frac_file='./input/GF3aTree.nc'		#frac of trees cover      - file path

   nl_tree_frac_file='./input/NTfrac_reorder_lat.nc'	#frac of needleleaf trees - file path
   tp_tree_frac_file='./input/tropfrac_reorder_lat.nc'	#frac of tropical   trees - file path
#  bl_tree_frac_file='./input/BLfrac_reorder_lat.nc'	#frac of broadleaf  trees - file_path

    ecotype_file='./input/EVT3b.nc'			#Gridded ecotypes ids     - file_path

     GtEcoEFfile="./db/GtEFbyEcotype.csv"		#Emission factor of each GT grouped by Ecotype
##Input files for BDPN: (not implemented yet)
        arid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_arid.nc'		# arid mask file
     nonarid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_non_arid.nc'     # non arid mask file
     landtype_files='./input/MEGAN31_Prep_Input_soil_191022/soil_landtype_'             # 24 landtype files
#      fert_files='input/MEGAN31_Prep_Input_soil_191022/soil_fert_'                     #
#nitro_depo_files='input/MEGAN31_Prep_Input_soil_191022/soil_nitrogen_mon'              #

#-------------------------------------------------
#(0) Get grid & proj parameters from GRIDDESC:
read projName xorig yorig dx dy nx ny nz <<< $( sed -n "/${GRIDNAME}/{n;p;q}" "$GRIDDESC_file" )
read COORDTYPE P_ALP P_BET P_GAM XCENT YCENT <<< $( sed -n "/${projName}/{n;p;q}" "$GRIDDESC_file" )
#read COORDTYPE truelat1 truelat2 stand_lon ref_lon ref_lat	#(en el namelist de wrf)
#truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
  if [ $COORDTYPE == 1 ]; then           #Geographic:
   srsOut="+proj=latlong +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 2 ]; then     #Lambert Conformal Conic:
   srsOut="+proj=lcc +lat_1=$P_ALP +lat_2=$P_BET +lon_0=$P_GAM +lat_0=$YCENT +a=6370000.0 +b=6370000.0 +units=m"
elif [ $COORDTYPE == 3 ]; then     #General Mercator
   srsOut="+proj=merc +lat_ts=$P_ALP +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 4 ]; then     #General tangent Stereografic
   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
elif [ $COORDTYPE == 5 ]; then     #UTM
   echo  "proyección: 5 (Universal Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 6 ]; then     #Polar Secant Stereographic
   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
elif [ $COORDTYPE == 7 ]; then     #Equatorial Mercator
   srsOut="+proj=merc +lat_ts=$P_ALP +a=6370000.0 +b=6370000.0"
elif [ $COORDTYPE == 8 ]; then     #Transverse Mercator
   echo  "proyección: 8 (Transverse Mercator) no soportada en esta aplicación."; stop
elif [ $COORDTYPE == 9 ]; then     #Lambert Azimuthal Equal-Area
   echo  "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada en esta aplicación."; stop
else
   echo  "codigo de proyección invalido. COORDTYPE"; stop
fi;

echo "SRS of Output Grids: $srsOut"

#xorig=-1970000;yorig=-2370000;
#srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 
#nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell

#Proyección: armar srsOut en base a los parametros de GRIDDESC
#srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 

#Grilla (sale de GRIDDESC)
xmin=$xorig; 
ymin=$yorig; 
#xmax=$(( ${xorig} *-1 )); 
#ymax=$(( ${yorig} *-1 ))
xmax=$( bc -l <<<"${xorig} *-1 " )
ymax=$( bc -l <<<"${yorig} *-1 " )
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

echo "Regridding input files for BDNP..."

echo " $arid_file     ->    arid.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF"  $arid_file     ./tmp_grids/arid.nc 
echo " $nonarid_file  -> nonarid.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF"  $nonarid_file  ./tmp_grids/nonarid.nc 

#Para Landtype hay 24 archivos(creo que son 24clases de land, y hay que ir sumandolos)
echo "   $landtype_files        -> landtype.nc   "
for LT in $(seq --format='%02.0f' 1 1 24)
do
	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF" ${landtype_files}${LT}* ./tmp_grids/landtype${LT}.nc 
done

#============================================
#Armo namelist input para prepmegan4cmaq.exe
cat << EOF > example.inp 
&control
start_date='${start_date}',
  end_date='${end_date}',

    griddesc_file='${GRIDDESC_file}'
        
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
!BDNP
         arid_file='./tmp_grids/arid.nc',     
        narid_file='./tmp_grids/nonarid.nc',  
           lt_file='./tmp_grids/landtype'  

/
EOF

#(2) Llamar a programa en Fortran que ensamble todas las grillas y hace los calculos pertinentes
#./prepmegan4cmaq.exe < example.inp






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
