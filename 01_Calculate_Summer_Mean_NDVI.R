#### CAFE Tutorial: Download and process NDVI
# Step 2: Calculate summer average

#Author: Allison James, Zach Popp

# Load packages
#
library(terra)
library(tidyverse)
library(sf)


# Set directory
#
datadir <- ".data/" #Fill in with your directory
interdir <- paste0(datadir, "masked/")
outdir <- paste0(datadir, "summer_mean/")

# List files (Aqua has much more data)
#
aqua_data_files <- list.files(path = datadir, pattern = "MYD13.*.hdf$")


# Read in data
#
data_summer <- rast(paste0(datadir, aqua_data_files))

# Assess characteristics
#
dim(data_summer)
res(data_summer) #231.6564 231.6564 about 250m
nlyr(data_summer)
names(data_summer)


########################### Quality control ###################################

r <- data_summer

# Explore the different layers in the raster (also done above). 
# We can see that layers 1 and 2 are the NDVI and EVI values, and 
# layer 3 contains the vegetation index "quality" layer.
names(r)

# For each tile (16-day snapshot), we need to mask areas that do not have
# sufficient pixel quality. For this step, we will break up our raster layers by
# time point and then mask the NDVI and EVI layers using the VI quality layer. 

nlayers_per_snap <- 12
nsnaps <- length(names(r)) / nlayers_per_snap
if (nsnaps != floor(nsnaps)) stop("Unexpected number of layers / snapshot size")

# The links below have more information about the different layers. 
# The VI quality layer has a binary value with several binary bits - each bit
# represents a different quality metric. 

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

# For more information, you can read:
# https://www.ctahr.hawaii.edu/grem/mod13ug/sect0005.html 
# https://lpdaac.usgs.gov/documents/621/MOD13_User_Guide_V61.pdf
# Page 18 has the description of QA data set.


# What we will do is make a matrix of certain bit ranges and values for those
# ranges that we want to mask out. Then we can use the modis_mask function 
# which takes the VI quality layer and this matrix as an input to mask the 
# NDVI and EVI layers. 

from <- c(1,9,11,12,16)
to   <- c(2,9,11,14,16)
reject <- c("01,10,11", "1", "1", "000,011,100,101,110,111", "1")
qa_bits <- cbind(from, to, reject)
qa_bits


for (i in seq_len(nsnaps)) {
  start <- (i-1)*nlayers_per_snap + 1
  end   <- start + nlayers_per_snap - 1
  message("Processing snapshot ", i, " (layers ", start, ":", end, ")")
  
  # Subset snapshot and simplify names
  r_snap <- r[[start:end]]
  names(r_snap) <- gsub('^\"|\"$', '', names(r_snap)) # remove surrounding quotes
  
  # Find QA layer (VI Quality)
  qc_idx <- grep("VI Quality", names(r_snap), ignore.case = TRUE)
  if (length(qc_idx) < 1) stop("Quality band not found in snapshot ", i)
  qc <- r_snap[[qc_idx[1]]]
  plot(qc, main = paste0("Quality snapshot ", i))
  
  # Build quality mask (16-bit QA)
  quality_mask <- modis_mask(qc, 16, qa_bits)
  plot(quality_mask, main = paste0("Quality mask snapshot ", i))
  
  # Choose which bands to mask: NDVI and EVI are the first two
  r_to_mask <- r_snap[[1:2]]
  
  # If modis_mask returns TRUE for good pixels (keep), use directly.
  # If it marks bad pixels as TRUE, uncomment the invert line below:
  # quality_mask <- !quality_mask
  
  # apply mask (mask() will set values to NA where quality_mask is NA/0)
  rmasked <- mask(r_to_mask, quality_mask)
  
  # write snapshot out as GeoTIFF (one file per snapshot)
  out_fname <- file.path(paste0(interdir, "summer_snapshot_", i, "_masked.tif"))
  writeRaster(rmasked, out_fname, overwrite = TRUE)
  message("Wrote: ", out_fname)

}


####################### Processing to Season ###################################

data_masked_files <- list.files(path = interdir,
                                pattern = "summer_snapshot*")

data_masked <- rast(paste0(interdir, data_masked_files))

# Now we have 6 layers: NDVI and EVI for each snapshot
names(data_masked)

# Add time to data. We need the original file names for this since they have
# the date in them.
#
stopifnot(unique(substr(varnames(data_summer), 10, 16)) != length(data_masked)/2)

# Save years and Julian days as a lists of numbers
yjjj_all <- substr(varnames(data_summer), 10, 16)
yjjj_unique <- unique(yjjj_all)

years <- as.integer(substr(yjjj_unique, 1, 4))
julians <- as.integer(substr(yjjj_unique, 5, 7))

# Convert to Date: Julian day 1 => Jan 1, so subtract 1 for origin offset
dates_unique <- mapply(function(y, j) as.Date(j - 1, origin = paste0(y, "-01-01")),
                       years, julians, SIMPLIFY = TRUE)

# Coerce to date format
data_dates <- as.Date(rep(dates_unique, each = 2), origin = "1970-01-01")

data_dates

# Assign the date
time(data_masked) <- data_dates


# Get EVI and NDVI layers only
#
ndvi_summer <- subset(data_masked, names(data_masked) == 
                             "250m 16 days NDVI")
evi_summer <- subset(data_masked, names(data_masked) == 
                            "250m 16 days EVI" )

# Add date to layer names so they are unique
#
current_names_ndvi <- gsub('^\"|\"$', '', names(ndvi_summer))
new_names <- paste0(current_names_ndvi, "_", format(unique(data_dates), "%Y%m%d"))
names(ndvi_summer) <- new_names
names(ndvi_summer)

current_names_evi <- gsub('^\"|\"$', '', names(evi_summer))
new_names <- paste0(current_names_evi, "_", format(unique(data_dates), "%Y%m%d"))
names(evi_summer) <- new_names
names(evi_summer)


# Bring in your administrative boundaries again.
# Example: Suffolk County in Massachusetts.
# This time we will use it to crop our data down to the boundaries. Right now we 
# have data for the full tile(s) that intersect with our area of interest, so we can use
# the admin boundaries to remove the excess pixels in the raster layers.
# Shapefiles for USA geometries can be found at: 
# https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
#

#ma_county <- st_read('./tl_us_county.shp') #Shapefile of all US Counties
#ma_county_suff <- ma_county[ma_county$GEOID == "25025", ] #Suffolk County
# ma_county_suff <- st_read(paste0(out_dir, "input_data/", "suffolk.shp")) 
# ma_county_suff <- vect(ma_county_suff)


# Change for your file
admin_boundaries <- st_read(paste0(in_dir, "name_of_your_shapefile.shp")) 

admin_boundaries <- vect(admin_boundaries)


# Project boundaries to crs of NDVI data
#
admin_boundaries_proj <- project(admin_boundaries, crs(data_masked[[1]]))


# Crop NDVI and EVI to administrative boundaries
#
crop_ndvi_summer <- crop(ndvi_summer, ext(admin_boundaries_proj))
crop_evi_summer <- crop(evi_summer, ext(admin_boundaries_proj))


# Plot this new cropped layer.
# We should now see that we no longer see the entire tile.
plot(crop_ndvi_summer[[1]])


# Calculate the annual summertime mean
#

annual_mean_by_year <- function(rast, years){
  
  nm <- names(rast)
  date_strs <- sub(".*(\d{8})$", "\1", nm)
  dates <- as.Date(date_strs, "%Y%m%d")
  out_list <- list()
  
  for (yr in years) {
    sel <- which(format(dates, "%Y") == as.character(yr))
    if (length(sel) == 0) {
      warning("No layers found for year ", yr)
      next
    }
    
    # subset and compute mean
    subr <- rast[[sel]]
    yr_mean <- app(subr, mean, na.rm = TRUE)
    names(yr_mean) <- paste0(names(rast)[1] %>% sub("\d{8}$", "", .), "_", yr) # base name + year
    out_list[[as.character(yr)]] <- yr_mean
  }
}

years <- 2020:2022  #Change for your years of interest
crop_ndvi_avgsummer_byyear <- annual_mean_by_year(crop_ndvi_summer, years)
crop_evi_avgsummer_byyear  <- annual_mean_by_year(crop_evi_summer, years)

plot(crop_ndvi_avgsummer_byyear)
plot(crop_evi_avgsummer_byyear)


writeRaster(crop_ndvi_avgsummer_byyear, paste0(outdir, "suffolk_mean_ndvi.tif"), 
            overwrite = T)
writeRaster(crop_evi_avgsummer_byyear, paste0(outdir, "suffolk_mean_evi.tif"), 
            overwrite = T)



