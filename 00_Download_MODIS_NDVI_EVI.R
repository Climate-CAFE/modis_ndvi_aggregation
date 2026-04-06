#### CAFE Tutorial: Download and process NDVI
#Author: Allison James, Zach Popp

# Adapted from this tutorial: https://rspatial.org/modis/2-download.html, which
# shares information on how to apply methods found in the book 
# “Geographic Information Analysis” by David O’Sullivan and David J. Unwin (2nd Edition, 2010).


############ Step 1: Download MODIS data from NASA EarthData ###############


# Load packages
#
library(terra)
library(tidyverse)
library(sf)
#remotes::install_github("rspatial/luna") #Use this line of code to install
# the luna package for the first time.
library(luna)
library(tigris) # Only necessary if you intend to download US census boundaries.

############ Set up ###########

# Downloading MODIS data from NASA Earthdata requires making an account. 
# You can do so here: https://urs.earthdata.nasa.gov/
# Your login credentials will need to be referenced when downloading the data,
# so it is convenient to store your username and password in text files. 

# Read in username and password
#
user <- readLines("./user.txt")
# Ensure password is saved in a secure location
pass <- readLines("./pass.txt")

# Set output directory for raster files
#
in_dir <- "./input_data/" # Where administrative boundaries are saved
out_dir <- "./modis_layers/"

# Bring in data for your administrative boundaries as a shapefile. 
# As an example, we've included Suffolk County in Massachusetts. 
# We will use the overall county boundary, then aggregate at the tract level.
# Shapefiles for USA geometries can be found at: 
# https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
#
# us_counties <- counties(year = 2020)
# st_write(us_counties, paste0(out_dir, "input_data/", "counties.shp")) # Save for easy loading

# us_counties <- st_read(paste0(out_dir, "input_data/", "counties.shp")) # Read back in

# ma_county_suff <- ma_county[ma_county$GEOID == "25025", ] # Suffolk County only
# st_write(ma_county_suff, paste0(out_dir, "input_data/", "suffolk.shp")) # Save for easy loading
# ma_county_suff <- vect(ma_county_suff)

# Replace file name your shapefile name:
admin_boundaries <- st_read(paste0(in_dir, "name_of_your_shapefile.shp")) 

admin_boundaries <- vect(admin_boundaries)

# Alternatively, you could use a bounding box instead:
# bounding_box <- data.frame(x = c(), y = c()) #Create df of coordinates
# admin_boundaries <- st_as_sf(bounding_box, coords = c("x", "y"), crs = 4326)


# Select a product- we have chosen "MODIS/Aqua Vegetation Indices 
# 16-Day L3 Global 250m SIN Grid V061". Read more below:
#https://www.earthdata.nasa.gov/data/catalog/lpcloud-myd13q1-061

product_aqua <- "MYD13Q1" # MODIS Aqua. 16-day, 250 m resolution; NDVI/EVI

# Set time frame
#
years <- 2020:2022 #Change for your years of interest


############ Download NDVI for each year ###########

# Use info to download data from NASA (Aqua).
# There will be one file for each 16-day period between the start and end date


for(yr in years){
  
  # Default option here is to download June 1 - Aug 31 for each year. 
  #If you have a specific date range, modify these lines.
  start <- paste0(yr, "-06-01") 
  end <- paste0(yr, "-08-31")
  
  tryCatch({
    luna::getNASA(product = product_aqua,
                  start = start,
                  end = end,
                  aoi = ext(admin_boundaries),
                  download = TRUE,
                  overwrite = TRUE,
                  path = paste0(out_dir, "modis/"),
                  username = user, #Reference user text file
                  password = pass) #Reference password text file
    message(glue::glue("Completed download for {yr}."))
  }, error = function(e) {
    message(glue::glue("Download failed for {yr}"))
  })
  
}

# The file naming convention of the downloaded files is as follows:
# Starts with the product name (MYD13Q1)
# Followed by the Julian date of acquisition formatted as AYYYYDDD, tile identifier,
# version of data collected (061), Julian date and time as production formatted as 
# YYYYDDDHHMMSS, and the data format (hdf).

# One of the file names for Suffolk County looks like this:
# "MYD13Q1.A2020169.h12v04.061.2020340091454.hdf"

# Bring in one of the rasters (replace with one of your files):
#
r_in <- rast(paste0(out_dir, "modis/", "one_of_the_files.hdf"))

# Plot it
#
# Disappointing? apply some stretching; see `?plotRGB` for more options
plotRGB(r_in, r = 1, g = 4, b = 3, stretch="lin")

