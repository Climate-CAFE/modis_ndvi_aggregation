# Created by: Zach Popp
# Adapted by: Allison James
# Date Created: 03/02/2025
# Version Number: v1
# Date Modified: 
# Modifications:
#   Code modified from: https://github.com/Climate-CAFE/modis5-daily-heat-aggregation
#
# *************************************************************** #
# ~~~~~~~  MODIS Re-Analysis Raster Aggregation Point      ~~~~~~~ #
# *************************************************************** #
#   Adapted from scripts developed by Keith Spangler, Muskaan Khemani for 
#       processing raster data onto polygon boundaries:
#       https://github.com/Climate-CAFE/population_weighting_raster_data/blob/main/Population_Weighting_Raster_to_Census_Units.R
#
# Overview:
#     Process MODIS rasters to administrative boundaries. This script is the 
#     second in a two-step raster processing process (following the steps to 
#     create mean summer NDVI and EVI in scripts 1-2). We will use the linked
#     grid cell - polygon points to extract NDVI and EVI and
#     then estimate polygon-level area weighted summer measures.
#
# Install required packages
#
library("terra")  # For raster data
library("sf")     # For vector data
library("plyr")   # For data management
library("doBy")   # For aggregation of data across groups
library("tidyverse") # For data management
library("lwgeom")
library("weathermetrics")
library("lubridate")
library("lutz")

sf_use_s2(FALSE) 
# S2 is for computing distances, areas, etc. on a SPHERE (using
# geographic coordinates, i.e., lat/lon in decimal-degrees); no need for this
# extra computational processing time if using PROJECTED coordinates,
# since these are already mapped to a flat surface. Here, pm25
# is indeed in geographic coordinates, but the scale of areas we are 
# interested in is very small, and hence the error introduced by 
# ignoring the Earth's curvature over these tiny areas is negligible and
# a reasonable trade off given the dramatic reduction in processing time. Moreover,
# the areas we calculate are not an integral part of the process
# and any error in that step would not materially impact the final output

# Check package version numbers
#
if (packageVersion("terra") < "1.7.78"   | packageVersion("sf") < "1.0.7" | 
    packageVersion("plyr")  < "1.8.7"    |
    packageVersion("doBy")  < "4.6.19"   | packageVersion("lwgeom") < "0.2.8") {
  cat("WARNING: packages are outdated and may result in errors.") }

################### User Defined Parameters ###################################

# Set up directories to read in and output data
#
modis_interdir <- "/Users/alliej/Documents/cafe/modis_ndvi_data/Suffolk/summer_mean/" # Directory where NDVI and EVI are output
outdir <- "/Users/alliej/Documents/cafe/modis_ndvi_data/OutputData/" # Output directory for tract-aggregated ndvi/evi
points_dir <- "/Users/alliej/Documents/cafe/modis_ndvi_data/MODIS_Fishnet/" # Input directory for extraction points


################### Run Processing to Aggregate ##########################

# Set name of geographic feature across which you want to summarize. This
# should be a column in the extraction_pts datasets
#
agg_geo <- "GEOID"

# Set year to process
#
years_to_agg <- c(2020)

# Define a function to calculate sums such that if all values are NA then it returns
# NA rather than 0.
#
sumfun <- function(x) { return(ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))) }

# Read in extraction points
#

region_in <- "25025" # Change for your GEOID
extraction_pts <- vect(paste0(points_dir, "MODIS_county_usreg", region_in, "_extraction_pts.gpkg"))


# Set up path for all MODIS NDVI and EVI files
#
modis_files <- list.files(path = modis_interdir, pattern = ("ndvi.*.tif$|evi.*.tif$"))

ndvi_all <- NULL
evi_all <- NULL
  
for (year in c(years_to_agg)) {
  
  cat("Now processing ", year, "\n")
  
  # Subset to year from all files (if applicable)
  #
  modis_files_yr <- modis_files[grepl(year, modis_files)]
  
  
  # Read in NDVI and EVI
  #
  ndvi <- rast(paste0(modis_interdir, modis_files_yr[grepl("ndvi", modis_files_yr)]))
  evi <- rast(paste0(modis_interdir,modis_files_yr[grepl("evi", modis_files_yr)]))
  
  # Rename the layers from "mean" to something more descriptive
  #
  stopifnot(length(names(ndvi)) == 1)
  stopifnot(length(names(evi)) == 1)
  names(ndvi) <- paste0("summer_mean_ndvi_", year)
  names(evi) <- paste0("summer_mean_evi_", year)
  
  # Reproject the layers to WGS84
  #
  ndvi <- terra::project(ndvi, "EPSG:4326")
  evi <- terra::project(evi, "EPSG:4326")
  
  if (is.null(ndvi_all)) {
    ndvi_all <- ndvi
  } else {
    ndvi_all <- c(ndvi_all, ndvi)   # combine layers (time dimension)
  }
  if (is.null(evi_all)) {
    evi_all <- evi
  } else {
    evi_all <- c(evi_all, evi)
  }
  
  list_rasters <- c(ndvi_all, evi_all)
  
  
    # Project points to WGS84 (coordinate system for ERA stack)
    #
    extraction_pts_cut <- project(extraction_pts, crs(list_rasters[[i]]))
    
    # Extract daily summaries to point-based grid
    #
    ndvi_pts <- terra::extract(list_rasters[[1]], extraction_pts_cut)
    evi_pts <- terra::extract(list_rasters[[2]], extraction_pts_cut)
    
    # Join results with extraction points (includes geographic identifiers)
    #
    ndvi_pts <- cbind(extraction_pts_cut, ndvi_pts)
    evi_pts <- cbind(extraction_pts_cut, evi_pts)
    
    ####################### Link Missing Data ##################################
    # Note that here were are linking points that are missing to the nearest
    # value. This may not be advisable in all contexts! Exploring where missing
    # points are is important to determine whether leaving as missing or joining
    # to nearest available points is advisable
    #
    # Use package sf to join points on coast and not covered by ERA5 land to
    # nearby temperature measures
    #
    ndvi_pts_sf <- st_as_sf(ndvi_pts)
    evi_pts_sf <- st_as_sf(evi_pts)
    
    # Examine where NA values are - for Suffok County, they are mostly water 
    # or very close to the coastline.
    #
    ndvi_pts_plot  <- ndvi_pts_sf %>%
      mutate(ndvi_status = if_else(is.na(summer_mean_ndvi_2020), "NA", "Not NA")) %>%
      ggplot() +
      geom_sf(aes(color = ndvi_status)) +
      scale_color_manual(values = c("NA" = "red", "Not NA" = "blue"))
    
    ndvi_pts_plot
    
    # Identify a date in year to assess missingness for points outside
    # land extent
    #
    varname_ndvi <- paste0("summer_mean_ndvi_", year)
    varname_evi<- paste0("summer_mean_evi_", year)
    
    # Subset missing data
    #
    ndvi_pts_missing <- ndvi_pts_sf[is.na(ndvi_pts_sf[[varname_ndvi]]), ]
    evi_pts_missing <- evi_pts_sf[is.na(evi_pts_sf[[varname_evi]]), ]
    
    # Subset available data
    #
    ndvi_pts_avail <- ndvi_pts_sf[!is.na(ndvi_pts_sf[[varname_ndvi]]), ]
    evi_pts_avail <- evi_pts_sf[!is.na(evi_pts_sf[[varname_evi]]), ]
    
    # Rejoin updated data with those already including temp
    #
    ndvi_pts <- rbind(ndvi_pts_avail, ndvi_pts_missing)
    evi_pts <- rbind(evi_pts_avail, evi_pts_missing)
    
    # Convert the extracted data to a data frame
    #
    ndvi_pts_df <- as.data.frame(ndvi_pts)
    evi_pts_df <- as.data.frame(evi_pts)
    
    # Remove geometry column
    #
    ndvi_pts_df$geometry <- NULL
    evi_pts_df$geometry <- NULL
    
    # Extract columns relevant to MODIS data
    #
    ndvi_cols <- names(ndvi_pts_df)[!names(ndvi_pts_df) %in% c("ID.1", names(extraction_pts))]
    evi_cols <- names(evi_pts_df)[!names(evi_pts_df) %in% c("ID.1", names(extraction_pts))]
    
    # Set names for MODIS variables based on input raster naming
    #
    ndvi_name <- paste0(substr(names(list_rasters[[1]])[1], 1, 16))
    evi_name <- paste0(substr(names(list_rasters[[2]])[1], 1, 15))
    
    # Transpose the data frame to get time series format (long-form)
    #
    ndvi_long <- ndvi_pts_df %>%
      pivot_longer(cols = all_of(ndvi_cols), names_to = "date", values_to = ndvi_name) 
    
    ndvi_long <- ndvi_long %>% mutate(date = str_extract(date, "\\d{4}$"))
    
    evi_long <- evi_pts_df %>%
      pivot_longer(cols = all_of(evi_cols), names_to = "date", values_to = evi_name) %>%
      select(UniqueID, date, !!sym(evi_name))
    
   evi_long <- evi_long %>% mutate(date = str_extract(date, "\\d{4}$"))
    
    # Combine data into single dataframe
    #
    MODIS_complete <- left_join(ndvi_long, evi_long, by = c("UniqueID", "date"))
  

############################## County Aggregation ############################
# In this step, we will use the extraction points to extract the MODIS grid cell
# underlying each portion of a ward across the entire raster stack of values.
# We again will follow the same process for each individual variable, and use 
# a loop to conduct the processing
#
varnames <- names(MODIS_complete)[!names(MODIS_complete) %in% c(names(extraction_pts), "date", "ID.1")]

for (i in 1:length(varnames)) {
  
  cat("Now processing ", varnames[i], "\n")
  
  # Before we calculate the final weighted average of the MODIS measure, we need to check for missing data.
  # If a value is NA on one of the polygons, then it will give an underestimate
  # since the weights will no longer add to 1. Example: there are two
  # polygons, each with 50% area. If Tmax is 30 C in one and NA in the other, then
  # the area weighted average (which removes NA values) would give: (30 * 0.5) + (NA * 0.5) = 15 C.
  # Therefore, we need to re-weight the weights based on the availability of data.
  #
  eqn <- as.formula(paste0("SpatWt ~ ", agg_geo, " + date")) 
  avail <- summaryBy(eqn,
                     data = MODIS_complete[which( !(is.na(MODIS_complete[[varnames[i]]])) ),],
                     FUN = sumfun)
  
  # Merge this value back into the longform MODIS data
  #
  MODIS_complete <- merge(MODIS_complete, avail, by = c(agg_geo, "date"), all.x = TRUE)
  
  # Re-weight the area weight by dividing by total available weight
  #
  MODIS_complete$SpatWt <- MODIS_complete$SpatWt / MODIS_complete$SpatWt.sumfun
  
  # QC: check that the weights of *available data* all add to 1
  #
  eqn <- as.formula(paste0("SpatWt ~ ", agg_geo, " + date"))
  check <- summaryBy(eqn,
                     data = MODIS_complete[which( !(is.na(MODIS_complete[[varnames[i]]])) ),],
                     FUN = sumfun)
  
  if (length(which(round(check$SpatWt.sumfun, 4) != 1)) > 0) {
    cat("ERROR: weights do not sum to 1", "\n"); break 
  } else {
    cat(":) weights sum to 1", "\n")
    MODIS_complete$SpatWt.sumfun <- NULL
  }
  
  # Multiply the variable of interest (here "newvarname") by the weighting value and then
  # sum up the resultant values within admin boundaries. This is an area-weighted average.
  # We will also divide by 100 million, as NDVI/EVI values range from -1 to 1 but are stored
  # as large whole numbers to preserve decimal points.
  #
  greenvar <- paste0(varnames[i], "_Wt")
  MODIS_complete[[greenvar]] <- MODIS_complete[[varnames[i]]] * MODIS_complete[["SpatWt"]] /100000000
  
  eqn <- as.formula(paste0(greenvar, " ~ ", agg_geo, " + date")) 
  
  final <- summaryBy(eqn, data = MODIS_complete, FUN = sumfun)
  
  # Automated QC to confirm that the dimensions are correct
  #
  if ( length(unique(as.data.frame(extraction_pts)[[agg_geo]]))  * length(unique(MODIS_complete$date)) != dim(final)[1]) {
    cat("ERROR: incorrect dimensions of final df", "\n"); break
  } else { cat(":) dimensions of final df are as expected", "\n") }
  
  # Set name for output
  #
  names(final)[grep(paste0("^", varnames[i]), names(final))] <- varnames[i]
  
  if (i == 1) {
    finaloutput <- final
  } else {
    finaloutput <- left_join(finaloutput, final, by = c(agg_geo, "date"))
  }
  
}

####################### Quality Control of Output ############################
cat("The final output has", dim(finaloutput)[1], "rows. \n")
cat("The first few lines of the output are: \n")
print(head(finaloutput))

# Automated QC: missing data
#
missing_ndvi_mean <- which(is.na(finaloutput$summer_mean_ndvi))
missing_evi_mean <- which(is.na(finaloutput$summer_mean_evi))

if (length(missing_evi_mean) > 0 | length(missing_ndvi_mean) > 0) {
  cat("WARNING: Note the number of missing ward-summers by variable: \n")
  cat("EVI mean:", length(missing_evi_mean), "\n")
  cat("NDVI mean:", length(missing_ndvi_mean), "\n")
  
  cat("The first few lines of missing evi (if any) are printed below: \n")
  print(head(finaloutput[missing_evi_mean,]))
  
  cat("The first few lines of missing ndvi (if any) are printed below: \n")
  print(head(finaloutput[missing_ndvi_mean,]))
  
} else { cat(":) No missing greenness values! \n") }

# Automated QC: impossible greenness values
#
num_temp_errors_evi <- length(which(finaloutput$summer_mean_evi < -1 |
                                       finaloutput$summer_mean_evi > 1))

num_temp_errors_ndvi <- length(which(finaloutput$summer_mean_ndvi < -1 |
                                       finaloutput$summer_mean_ndvi > 1))

if (num_temp_errors_evi > 0 | num_temp_errors_ndvi > 0 ) { 
  
  print("ERROR: impossible greenness values. Applicable rows printed below:")
  print(finaloutput[which(finaloutput$summer_mean_evi < -1 |
                            finaloutput$summer_mean_evi > 1),])
  print(finaloutput[which(finaloutput$summer_mean_ndvi < -1 |
                            finaloutput$summer_mean_ndvi > 1),])
  
} else { print(":) all greenness values are of correct *relative* magnitude") }

# Output results by year to output directory
#
saveRDS(finaloutput, paste0(outdir, "/", "county_agg_MODIS_", year, "_ndvi_evi_region", region_in, ".rds"))
write_csv(finaloutput, paste0(outdir, "/", "county_agg_MODIS_", year, "_ndvi_evi_region", region_in, ".csv"))

}
