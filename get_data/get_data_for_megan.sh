#!/bin/bash

#=================================
#Leaf Area Index (LAI) 
#is used to quantify the amount of foliage in a landscape. MEGAN uses a value, referred to as LAIv),  that is representative of the m2 of one sided leaf area per m2 of canopy area (vegetation covered area).
# (1)  LAIv = LAI/VCF * 10
#the unit used in MEGAN3.2 pre-processor is m2 leaf area per 10 m2 surface area so if you generate LAIv in m2/m2 then multiply this by 10 to get m2 leaf area per 10 m2.
#Namming convention: laiv$YYYY$MM.nc

#Sources of data:
# Leaf Area Index (LAI)
#	- MODIS (MCD15A2.005)    https://e4ftl01.cr.usgs.gov/MOTA/MCD15A2H.006/
#	- Landsat
#	- Sentinel
#
# Vegetation cover fraction (VCF)
#	- MODIS	(MOD44B)	https://e4ftl01.cr.usgs.gov/MOLT/MOD44B.006/
#	- Landsat
#	- Sentinel 
#	- AVHRR


#They are based on MODIS (MOderate Resolution Imaging Spectroradiometer) satellite data using the LAI estimates of Zhang et al JGR (2004) doi:10.1029/2004JD004720 and the vegetation cover fraction of Hansen et al.  Earth Interactions 7: 1-15 (2003).

#They have a file naming convention of laiv2003MM where MM is month of the year. Urban and wetland areas and other regions with missing LAI  data in this MODIS database were assigned an LAIv based on the average for the surrounding region.

#Each 1-year dataset is 31MB zipped and are based on MODIS (MOderate Resolution Imaging Spectroradiometer) satellite version 5 product (MCD15A2.005) and maximum green vegetation fraction from USGS, which is also based on MODIS.  

#=================================
#Growth Forms

#Crops
#Shurbs
#Grass
#Tree	(needleleaf, broadleaf, tropical)

#Sources:
# + Guenther et al. 2012 (MEGAN2 viene con sus datasets)
# + PROBA-V  Copernicus Global Land Covert (tiene frac covers de Forest, Shurblands, Herbacious, Croplands (Hasta 2019)

#=================================
#Ecotypes

#Source:
# + Guenther et al. 2012 (MEGAN2 viene con sus datasets)

#Me gustaria entender un poco mejor esta variable.



