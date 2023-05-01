# MEGAN Emissions

> `prepmegan4cmaq` is a fortran program that grids all the data needed to run biogenic emissions inside CMAQ model.

## Dependencies:
 +  Fortran GNU compiler
 +  NetCDF library
 +  GDAL/OGR programs (`sudo apt install gdal-bin`)

## Get MEGAN input data:

All the data required to run this pre-procesor is freely available from [UCI BAI webpage](https://bai.ess.uci.edu/megan/data-and-code/) and has been produced by the team of Alex Guenther.

Data required:
+ Leaf Area Index / Vegetation Cover Fraction (LAIv)
+ Growth Form (fraction): crops, grass, shurbs, tree.
+ Canopy type (fraction): tropical trees, needleleaf trees.
+ Ecotype
+ Soil data (for BDSNP soil NO algorithm): Fertilizer, Land Fraction, Climate data (arid/non-arid), Nitrogen deposition.


## Build
Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your `nc-config --libdir` and `nc-config --includedir`.

`> make`

If the compilation is successful, the executable `prepmegan4cmaq.exe` should be created.

## Run

Edit the header of `prepmegan4cmaq.sh` script that contains the following variables:

```shell
start_date="2019-01-01" #"%Y-%m-%d %H"
  end_date="2019-01-01" #"%Y-%m-%d %H"

   srsInp="epsg:4326"   #Cambiar <solo si> los archivos NO vienen en latlon.

GRIDDESC_file="./GRIDDESC_example"
GRIDNAME="Argentina"

    crop_frac_file='./input/GF3aCrop.nc'             #frac of crop cover       - file path
   grass_frac_file='./input/GF3aGrass.nc'            #frac of grass cover      - file path
   shrub_frac_file='./input/GF3aShurb.nc'            #frac of shurbs cover     - file path
    tree_frac_file='./input/GF3aTree.nc'             #frac of trees cover      - file path
 nl_tree_frac_file='./input/NTfrac_reorder_lat.nc'   #frac of needleleaf trees - file path
 tp_tree_frac_file='./input/tropfrac_reorder_lat.nc' #frac of tropical   trees - file path
#bl_tree_frac_file='./input/tropfrac_reorder_lat.nc' #frac of broadleaf  trees - file_path
      ecotype_file='./input/EVT3b.nc'                #(not so clear what it represents) - file_path
        laiv_files='./input/laiv2003'                #(LAIv=LAI/VEGCOVER) file path just the name before month indicator!

```

Note that the variables must be adjusted to match the appropriate values for your system.

Then execute `prepmegan4cmaq.sh`:

`> ./prepmegan4cmaq.sh` 

This script will regrid all the input files using GDAL and then execute the `prepmegan4cmaq.exe` program.

Please feel free to contact the developer if you have any issues or suggestions.


## Planned future improvements:
 + [ ] BDNP Fert support
 + [ ] Roboust GRIDDESC reader
 + [ ] Write temporal-dependent data only for the month within startdate and enddate

