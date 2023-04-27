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
      laiv_files='./input/laiv2003'			#(LAIv=LAI/VEGCOVER) file path just the name before month indicator!

      cropf_file='./input/GF3aCrop.nc' 		        #frac of crop cover       - file path
     grassf_file='./input/GF3aGrass.nc'		        #frac of grass cover      - file path
     shrubf_file='./input/GF3aShrub.nc'		        #frac of shrubs cover     - file path
      treef_file='./input/GF3aTree.nc'		        #frac of trees cover      - file path

   nl_treef_file='./input/NTfrac_reorder_lat.nc'	#frac of needleleaf trees - file path
   tp_treef_file='./input/tropfrac_reorder_lat.nc'	#frac of tropical   trees - file path
#  bl_treef_file='./input/BLfrac_reorder_lat.nc'	#frac of broadleaf  trees - file_path

    ecotype_file='./input/EVT3b.nc'			#(not so clear what it represents) - file_path

#For BDPN:
       arid_file='input/MEGAN31_Prep_Input_soil_191022/soil_climate_arid.nc'
    nonarid_file='input/MEGAN31_Prep_Input_soil_191022/soil_climate_non_arid.nc'
   landtype_file='input/MEGAN31_Prep_Input_soil_191022/soil_landtype_'
      fert_files='input/MEGAN31_Prep_Input_soil_191022/soil_fert_'
nitro_depo_files='input/MEGAN31_Prep_Input_soil_191022/soil_nitrogen_mon'
#       lai_file='./input/lai*'
#   veg_cov_file='./input/veg_cov*3'
#    wrfout_file='./input/wrfout_$YYYY$MM$DD_00:00:00_d01.nc'
#-------------------------------------------------
#(0) Get grid & proj parameters from GRIDDESC:
read projName xorig yorig dx dy nx ny nz <<< $( sed -n "/${GRIDNAME}/{n;p;q}" "$GRIDDESC_file" )
read COORDTYPE P_ALP P_BET P_GAM XCENT YCENT <<< $( sed -n "/${projName}/{n;p;q}" "$GRIDDESC_file" )

xorig=-1970000;yorig=-2370000;
truelat1=-50;truelat2=-20;stand_lon=-65;ref_lon=-65;ref_lat=-35;
srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 
nx=197;ny=237;nz=1;dx=20000;dy=20000;     #ncols,nrows,xcel, ycell

#Proyección: armar srsOut en base a los parametros de GRIDDESC
#srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m" 

#Grilla (sale de GRIDDESC)
xmin=$xorig; ymin=$yorig; xmax=$(( $xorig*-1 )); ymax=$(( $yorig*-1 ))
#------------------------------------------------
#(1) Re-gridding: Agarrar los inputs y regrillarlos (e interpolar) segun GRIDDESC.
#interpolation methods: near (default), bilinear, cubic, cubicspline, lanczos, average, rms, mode,  max, min, med, Q1, Q3, sum
if [ ! -d "./tmp_grids" ]; then
    mkdir tmp_grids
fi

echo "Regridding input files..."

#echo "   $cropf_file    -> crop.nc   "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $cropf_file ./tmp_grids/crop.nc 
#echo "   $grassf_file   -> grass.nc  "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $grassf_file ./tmp_grids/grass.nc 
#echo "   $shrubf_file   -> shrub.nc  "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $shrubf_file ./tmp_grids/shrub.nc
#echo "   $treef_file    -> tree.nc   "; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $treef_file ./tmp_grids/tree.nc 
#echo "   $nl_treef_file -> nl_tree.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $nl_treef_file ./tmp_grids/nl_tree.nc 
#echo "   $tp_treef_file -> tp_tree.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $tp_treef_file ./tmp_grids/tp_tree.nc
echo "   $ecotype_file      -> ecotype.nc"; gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r mode -f "NetCDF"      $ecotype_file ./tmp_grids/ecotype.nc 

#echo "   $laiv_files        -> laiv.nc   "
#for MM in $(seq --format='%02.0f' 1 1 12)
#do
#	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" ${laiv_files}${MM}* ./tmp_grids/laiv${MM}.nc 
#done

#Armo namelist input para prepmegan4cmaq.exe
cat << EOF > example.inp 
&control
start_date='${start_date}',
  end_date='${end_date}',

    griddesc_file='${GRIDDESC_file}'

    cropf_file='./tmp_grids/crop.nc', 
   grassf_file='./tmp_grids/grass.nc',
   shrubf_file='./tmp_grids/shrub.nc',
    treef_file='./tmp_grids/tree.nc', 
 nl_treef_file='./tmp_grids/nl_tree.nc',		
 tp_treef_file='./tmp_grids/tp_tree.nc',		
 bl_treef_file='./tmp_grids/bl_tree.nc',		

      ecotype_file='./tmp_grids/ecotype.nc',				

        laiv_files='./tmp_grids/laiv',				
/
EOF

#(2) Llamar a programa en Fortran que ensamble todas las grillas y hace los calculos pertinentes
#./prepmegan4cmaq.exe < example.inp

