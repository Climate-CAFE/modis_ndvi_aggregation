#### CAFE Tutorial: Download and process NDVI
# Step 1: Download MODIS data from NASA EarthData

#Author: Allison James, Zach Popp
#Date: 1/22/2026


# Load packages
#
library(terra)
library(tidyverse)
library(sf)
#remotes::install_github("rspatial/luna") #Use this line of code to install
# the luna package for the first time.
library(luna)
library(tigris)
options(tigris_use_cache = FALSE)


############ Set up ###########

# Read in username and password
#
user <- readLines("/Users/alliej/Desktop/bu/CAFE/CAFE_MODIS/user.txt")
# Ensure password is saved in a secure location
pass <- readLines("/Users/alliej/Desktop/bu/CAFE/CAFE_MODIS/pass.txt")

# Set directory
#
datadir <- "/Users/alliej/Documents/cafe/modis_ndvi_aggregation/example_outputs/" #Fill in with your directory
out_dir <- "/Users/alliej/Documents/cafe/modis_ndvi_data/"


# https://www.earthdata.nasa.gov/data/catalog/lpcloud-mod13q1-061#toc-file-naming-convention
product_aqua <- "MYD13Q1" # MODIS Aqua. 16-day, 250 m; NDVI/EVI

# Set time frame
#
start <- "2020-07-01"
end <- "2020-07-31"


# Bring in data for Suffolk County in Massachusetts. 
# Shapefiles for USA geometries can be found at: 
# https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
#
#ma_county <- counties(year = 2020)
#st_write(ma_county, paste0(out_dir, "input_data/", "counties.shp"))

#ma_county <- st_read('./tl_2020_us_county.shp') #Shapefile of all US Counties
ma_county <- st_read(paste0(out_dir, "input_data/", "counties.shp")) #Shapefile of all US Counties

ma_county_suff <- ma_county[ma_county$GEOID == "25025", ] #Suffolk County
ma_county_suff <- vect(ma_county_suff)


# Use info to download data from NASA (Aqua)
#
data_ndvi_aqua <- luna::getNASA(product_aqua, 
                                start, 
                                end, 
                                aoi=ext(ma_county_suff), 
                                download=TRUE, 
                                overwrite=TRUE,
                                path=paste0(out_dir, "Suffolk/"), 
                                username=user, #Reference user text file
                                password=pass) #Reference password text file


# Bring in raster
#
r_in <- rast(paste0(out_dir, "Suffolk/", "MYD13Q1.A2020169.h12v04.061.2020340091454.hdf"))

# Plot it
#
# Disappointing? apply some stretching; see `?plotRGB` for more options
plotRGB(r_in, r = 1, g = 4, b = 3, stretch="lin")

