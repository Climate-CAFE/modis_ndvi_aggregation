#### CAFE Tutorial: Download and process NDVI
# Step 2: Calculate summer average and aggregate

#Author: Allison James, Zach Popp
#Date: 2/4/2026


# Load packages
#
library(terra)
library(tidyverse)
library(sf)


# Set directory
#
datadir <- ".data/" #Fill in with your directory


# List files (Aqua has much more data)
#
aqua_data_files <- list.files(path = paste0(datadir, "Suffolk/"),
                              pattern = "MYD13.*.hdf$")



# Data are available in 16 day increments. Summer is most of interest. 
# Adjust dates for your time period. 
# If you are working with a period of multiple years, that includes a leap year,
# you may need to repeat this process for that year specifically, as summer
# will be a different date range.
#
format(as.Date("2020-06-01"), "%j")
format(as.Date("2020-08-31"), "%j")

aqua_summer_2020 <- aqua_data_files[substr(aqua_data_files, 14, 16) %in% c(153:244)]

# Read in 2018 data
#
data_summer_2020 <- rast(paste0(datadir, "Suffolk/", aqua_summer_2020))

# Assess characteristics
#
dim(data_summer_2020)
res(data_summer_2020) #231.6564 231.6564 about 250m
nlyr(data_summer_2020)
names(data_summer_2020) # 3*12 layers = 36
# 16 days: ndvi, evi, vi quality, red reflectance, nir reflec, blue reflec, 
# mir reflec, view zenith, sun zenity, relative azimuth, composite day, pixel reliablity


########################### Quality control ###################################


r <- data_summer_2020

#https://www.ctahr.hawaii.edu/grem/mod13ug/sect0005.html 
#https://lpdaac.usgs.gov/documents/621/MOD13_User_Guide_V61.pdf
#page 18 for the description of QA data set

from <- c(1,9,11,12,16)
to   <- c(2,9,11,14,16)
reject <- c("01,10,11", "1", "1", "000,011,100,101,110,111", "1")
qa_bits <- cbind(from, to, reject)
qa_bits

# Bits 1-2: Overall quality. Keep 00 (good) and maybe 01 (check other QA)
# Bits 3-6 (usefulness): Has scale of quality.
# Bits 7-8: Aerosol Quantity (climatology, low, intermediate, high)
# Bit 9: Adjacent clouds detected
# Bit 10: Atmosphere BRDF correction
# Bit 11: Mixed clouds
# Bits 12-14: Land/Water mask (ocean, land, coastline, shallow inland water, 
# ephemeral water, deep inland water, continental ocean, deep ocean). Probably want to 
# keep land 01 and 101 coastline.
# Bit 15: possible snow/ice (hopefully not relevant for summer)
# Bit 16: possible shadow

# Explore the different layers in the raster (also done above). 
# We can see that layers 1 and 2 are the NDVI and EVI values, and 
# layer 3 contains the vegetation index "quality" layer.
names(r)

qc <- r[[3]] #Band 3 has the VI Quality layer
plot(qc, main = "Quality")


# This MODIS mask function in the luna package can create a "mask" for each 
# layer based on our selection of bad quality pixels. It takes the number of 
# bits as an input so we should use the VI Quality layer. 
# The other arguments are a raster and a 3-column matrix. 
# The first two columns represent the start and end of a range of bits, and the
# third column contains the value(s) that should be rejected.
quality_mask <- modis_mask(qc, 16, qa_bits)
plot(quality_mask, main = "Quality mask")

rmask <- mask(r, quality_mask)

#### add a legend to mine ####
plotRGB(rmask, r = 2, g = 1, b = 4, main='False color composite', stretch="lin")


writeRaster(rmask, paste0(datadir, "Suffolk/", "suffolk_july_2020_masked.tif"), 
            overwrite = T)

####################### Processing to Season ###################################

data_masked <- rast(paste0(datadir, "Suffolk/", "suffolk_july_2020_masked.tif"))

# Add time to data
#
time(data_masked) <- as.numeric(substr(varnames(data_summer_2020), 10, 16))

# Convert Julian day to Date, note subtract 1 because origin is day 1 (Jan 1)
time(data_masked) <- as.Date(as.numeric(substr(varnames(data_summer_2020), 14, 16)) - 1, 
                             origin = paste0(substr(varnames(data_summer_2020), 10, 13), "-01-01"))

# Get EVI and NDVI layers only
#
ndvi_summer_2020 <- subset(data_masked, names(data_masked) == 
                             "\"250m 16 days NDVI\"")
evi_summer_2020 <- subset(data_masked, names(data_masked) == 
                            "\"250m 16 days EVI\"" )


# Bring in data for Suffolk County in Massachusetts again.
# This time we will use it to crop our data down to Suffolk County. Right now we 
# have data for the tile(s) that intersect with our area of interest, so we can use
# the outline of the county to remove the excess pixels in the raster layers.
# Shapefiles for USA geometries can be found at: 
# https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
#
ma_county <- st_read('./tl_2020_us_county.shp') #Shapefile of all US Counties
ma_county_suff <- ma_county[ma_county$GEOID == "25025", ] #Suffolk County
ma_county_suff <- vect(ma_county_suff)


# Try to project MA counties to crs of data
#
ma_county_suff_proj <- project(ma_county_suff, crs(data_masked[[1]]))
ext(ma_county_suff_proj)

# Crop NDVI to suffolk county
#
suff_ndvi_summer_2020 <- crop(ndvi_summer_2020, ext(ma_county_suff_proj))
suff_evi_summer_2020 <- crop(evi_summer_2020, ext(ma_county_suff_proj))


# Plot this new Suffolk county only layer.
# We should now see that we no longer see the entire tile, just Suffolk County.
plot(suff_ndvi_summer_2020[[1]])


# Calculate mean for summertime
#
suff_ndvi_avgsummer_20 <- app(suff_ndvi_summer_2020, mean)
suff_evi_avgsummer_20 <- app(suff_evi_summer_2020, mean)


writeRaster(suff_ndvi_avgsummer_20, paste0(datadir, "Suffolk/summer_mean/", "suffolk_2020_mean_ndvi.tif"), 
            overwrite = T)
writeRaster(suff_evi_avgsummer_20, paste0(datadir, "Suffolk/summer_mean/", "suffolk_2020_mean_evi.tif"), 
            overwrite = T)



