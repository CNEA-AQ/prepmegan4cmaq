# MEGAN Emissions

> `prepmegan4cmaq` is a fortran program that grids all the data needed to run biogenic emissions inside CMAQ model.

## Dependencies:
 +  Fortran GNU compiler
 +  NetCDF library

## Get MEGAN input data:

All the data required to run this pre-processor is freely available from [UCI BAI webpage](https://bai.ess.uci.edu/megan/data-and-code/) and has been produced by the team of Alex Guenther.

Data required:
+ Leaf Area Index / Vegetation Cover Fraction (LAIv).
+ Growth type (fraction): crop, grass, shurb, tree, tropical trees, needleleaf trees..
+ Ecotype.
+ Soil data (for BDSNP soil NO algorithm): Land type, Climate data (arid/non-arid), Nitrogen deposition and soil Nitrogen from fertilizers.

## Build
Go to the ``src`` directory:

`> cd src/`

Edit the Makefile to set the compiler and path to NetCDF lib and include files. Check your `nc-config --libdir` and `nc-config --includedir`.

`> make`

If the compilation is successful, the executable `prepmegan4cmaq.exe` should be created.

## Run

Edit the namelist `example.inp` that contains the following variables:

```fortran
&control
griddesc_file='./GRIDDESC'
     gridname='LCC_TAN_TEST',

ecotypes_file='veg_Ecotypes.nc',
growtype_file='veg_GrowthFormFracions.nc',
    laiv_file='veg_LAIv.nc',
 climate_file='soil_climate.nc',
   ferti_file='soil_fert.nc',
landtype_file='soil_landtype.nc',
   nitro_file='soil_nitro.nc',

GtEcoEF_file="./db/GtEFbyEcotype.csv"      !Emission factor of each GT grouped by Ecotype
/
```

Note that the variables must be adjusted to match the appropriate values for your system.

Then execute `prepmegan4cmaq.exe`:

`> ./prepmegan4cmaq.exe < example.inp` 

Please feel free to contact the developer if you have any issues or suggestions.

## Planned future improvements:

 + [x] Reduce number input files and output files. 
   - [x] Group LAI files into one.
   - [x] Group GT frac files into one.
   - [x] Group Fert files and Nitro files.
   - [x] Group landtype files by using int as id. 
 + [x] Input and output files and variables with minningfull names.
 + [x] Remove GDAL/OGR dependence 
   - [x] Roboust GRIDDESC reader.
   - [x] coordinate transformations functions
   - [x] interpolation subroutines (bilinear, bicubic, average, mode, median)
 + [ ] BDNP (Fert variable) support.
 + [ ] Add some scripts to download and prepare some input files (LAI, Fert, Nitro, etc.)

