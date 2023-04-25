#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#------------------------------------------------
#Input-Data:
start_date="2019-01-01"	#"%Y-%m-%d %H"
  end_date="2019-01-01"	#"%Y-%m-%d %H"

   srsInp="epsg:4326"	#Cambiar <solo si> los archivos NO vienen en latlon.

GRIDDESC_file="./GRIDDESC_example"
GRIDNAME="Argentina"

#Input Files:
    crop_frac_file='./input/GF3aCrop.nc' 			#frac of crop cover       - file path
   grass_frac_file='./input/GF3aGrass.nc'			#frac of grass cover      - file path
   shrub_frac_file='./input/GF3aShrub.nc'			#frac of shrubs cover     - file path
    tree_frac_file='./input/GF3aTree.nc'			#frac of trees cover      - file path
 nl_tree_frac_file='./input/NTfrac_reorder_lat.nc'		#frac of needleleaf trees - file path
 tp_tree_frac_file='./input/tropfrac_reorder_lat.nc'		#frac of tropical   trees - file path
# bl_tree_frac_file='./input/tropfrac_reorder_lat.nc'		#frac of broadleaf  trees - file_path
      ecotype_file='./input/EVT3b.nc'				#(not so clear what it represents) - file_path
        laiv_files='./input/laiv2003'				#(LAIv=LAI/VEGCOVER) file path just the name before month indicator!
#	  lai_file='./input/lai*'
#      veg_cov_file='./input/veg_cov*3'
#       wrfout_file='./input/wrfout_$YYYY$MM$DD_00:00:00_d01.nc'

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

#Date params:
#start_date_s=$(date -d "$start_date" +%s)
#  end_date_s=$(date -d "$end_date  " +%s)
#------------------------------------------------
#(1) Re-gridding: Agarrar los inputs y regrillarlos (e interpolar) segun GRIDDESC.
#interpolation methods: near (default), bilinear, cubic, cubicspline, lanczos, average, rms, mode,  max, min, med, Q1, Q3, sum
if [ ! -d "./tmp_grids" ]; then
    mkdir tmp_grids
fi

echo "Regridding input files..."

gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $crop_frac_file ./tmp_grids/crop.nc 
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $grass_frac_file ./tmp_grids/grass.nc 
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"   $shrub_frac_file ./tmp_grids/shrub.nc
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"    $tree_frac_file ./tmp_grids/tree.nc 
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $nl_tree_frac_file ./tmp_grids/nl_tree.nc 
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" $tp_tree_frac_file ./tmp_grids/tp_tree.nc
gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut" -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF"      $ecotype_file ./tmp_grids/ecotye.nc 

for MM in $(seq --format='%02.0f' 1 1 12)
do
	gdalwarp -q -overwrite -s_srs "$srsInp" -t_srs "$srsOut"  -te $xmin $ymin $xmax $ymax -tr $dx $dy -r bilinear -f "NetCDF" ${laiv_files}${MM}* ./tmp_grids/laiv${MM}.nc 
done

#Armo namelist input para prepmegan4cmaq.exe
cat << EOF > example.inp 
&control
start_date='${start_date}',
  end_date='${end_date}',

    griddesc_file='${GRIDDESC_file}'

    crop_frac_file='./tmp_grids/crop.nc', 
   grass_frac_file='./tmp_grids/grass.nc',
   shrub_frac_file='./tmp_grids/shrub.nc',
    tree_frac_file='./tmp_grids/tree.nc', 
 nl_tree_frac_file='./tmp_grids/nl_tree.nc',		
 tp_tree_frac_file='./tmp_grids/tp_tree.nc',		
 bl_tree_frac_file='./tmp_grids/bl_tree.nc',		

      ecotype_file='./tmp_grids/ecotype.nc',				

        laiv_files='./tmp_grids/laiv',				
/
EOF

#(2) Llamar a programa en Fortran que ensamble todas las grillas y hace los calculos pertinentes
#./prepmegan4cmaq.exe < example.inp

