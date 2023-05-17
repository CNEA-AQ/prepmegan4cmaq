# MEGAN Emissions

> `prepmegan4cmaq` is a fortran program that grids all the data needed to run biogenic emissions inside CMAQ model.

## Dependencies:
 +  Fortran GNU compiler
 +  NetCDF library
 +  GDAL/OGR programs (`sudo apt install gdal-bin`)

## Get MEGAN input data:

All the data required to run this pre-processor is freely available from [UCI BAI webpage](https://bai.ess.uci.edu/megan/data-and-code/) and has been produced by the team of Alex Guenther.

Data required:
+ Leaf Area Index / Vegetation Cover Fraction (LAIv)
+ Growth Form (fraction): crop, grass, shurb, tree.
+ Canopy type (fraction): tropical trees, needleleaf trees.
+ Ecotype
+ Soil data (for BDSNP soil NO algorithm): Land type, Climate data (arid/non-arid), Nitrogen deposition and soil Nitrogen from fertilizers.


## Build
Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your `nc-config --libdir` and `nc-config --includedir`.

`> make`

If the compilation is successful, the executable `prepmegan4cmaq.exe` should be created.

The command `gdalwarp` (from GDAL/OGR) is used to regrid, interpolate and reproject the inputs grids, you can get it with `sudo` in ubuntu linux distributions by:

`> sudo apt install gdal-bin`

## Run

Edit the header of `prepmegan4cmaq.sh` script that contains the following variables:

```shell
#Input-Data:
start_date="2019-01-01" #"%Y-%m-%d %H"
  end_date="2019-01-01" #"%Y-%m-%d %H"

  srsInp="epsg:4326"    #(latlon spatial reference system of input files)

GRIDDESC_file="/home/usuario/runs/papila2019/cmaq/mcip/GRIDDESC" #"./GRIDDESC_example"
GRIDNAME="PAPILAGRID" #"Argentina"

#Input Files:
       laiv_files='./input/laiv2003'                     #LAIv=LAI/VegCover  (-) -  netCDFfiles path PREFIX

   crop_frac_file='./input/GF3aCrop.nc'                  #frac of crop   cover (%) - netCDF file path
  grass_frac_file='./input/GF3aGrass.nc'                 #frac of grass  cover (%) - netCDF file path
  shrub_frac_file='./input/GF3aShrub.nc'                 #frac of shrubs cover (%) - netCDF file path
   tree_frac_file='./input/GF3aTree.nc'                  #frac of trees  cover (%) - netCDF file path

nl_tree_frac_file='./input/NTfrac_reorder_lat.nc'        #frac of needleleaf trees (%) - netCDF file path
tp_tree_frac_file='./input/tropfrac_reorder_lat.nc'      #frac of tropical   trees (%) - netCDF file path

     ecotype_file='./input/EVT3b.nc'                     #Gridded ecotypes ids         - netCDF file_path

     GtEcoEFfile="./db/GtEFbyEcotype.csv"                #Emission factor of each GT grouped by Ecotype

#Input files for BDSNP:
arid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_arid.nc'          #     arid mask (0/1)- netCDF file path
nonarid_file='./input/MEGAN31_Prep_Input_soil_191022/soil_climate_non_arid.nc'   # non arid mask (0/1)- netCDF file path
  landtype_files='./input/MEGAN31_Prep_Input_soil_191022/soil_landtype_'         # landtype           - netCDF files PREFIX

  nitro_depo_files='input/MEGAN31_Prep_Input_soil_191022/soil_nitrogen_mon'      # soil-NO deposition of each month (kg/m2/s) - netCDF files PREFIX
  #     fert_files='input/MEGAN31_Prep_Input_soil_191022/soil_fert_'             #Reservoir of N associated w/ manure and fertilizer (mg/m3) - netCDF files PREFIX
```

Note that the variables must be adjusted to match the appropriate values for your system.

Then execute `prepmegan4cmaq.sh`:

`> ./prepmegan4cmaq.sh` 

This script will regrid all the input files using GDAL and then execute the `prepmegan4cmaq.exe` program.

Please feel free to contact the developer if you have any issues or suggestions.


## Planned future improvements:

 + [ ] BDNP (Fert variable) support.
 + [ ] Reduce number input files and output files. 
   - [ ] Group LAI files into one.
   - [ ] Group GT frac files into one.
   - [ ] Group Fert files and Nitro files.
   - [ ] Group landtype files by using int as id. 
 + [ ] Input and output files and variables with minningfull names.
 + [ ] Remove GDAL/OGR dependence 
   - [x] Roboust GRIDDESC reader.
   - [x] coordinate transformations functions
   - [ ] interpolation subroutines
 + [ ] Add some scripts to download and prepare some input files (LAI, Fert, Nitro, etc.)

