# Created by: Zach Popp
# Date Created: 05/15/2025
# Version Number: v1
# Date Modified: 
# Modifications:
#   Code modified from: https://github.com/Climate-CAFE/era5-daily-heat-aggregation
#
# *************************************************************** #
# ~~~~~~~  ERA5 Re-Analysis Raster Aggregation Point      ~~~~~~~ #
# *************************************************************** #
#   Adapted from scripts developed by Keith Spangler, Muskaan Khemani for 
#       processing raster data onto polygon boundaries:
#       https://github.com/Climate-CAFE/population_weighting_raster_data/blob/main/Population_Weighting_Raster_to_Census_Units.R
#
# Overview:
#     Process ERA5 rasters to administrative boundaries. This script is the 
#     second in a two-step raster processing process (following the steps to 
#     create hourly WBGT and UTCI in scripts 1-4). We will use the linked
#     grid cell - polygon points to extract WBGT and UTCI daily values and
#     then estimate polygon-level area weighted daily measures.
#
#     First, we need to estimate daily minimum, mean, and maximum WBGT and UTCI
#     from the hourly data. To get this we need to convert from the UTC time zone
#     in which the ERA5 data is distributed, to local time (as the min, mean, 
#     and max should be based on local time). We are here applying the code
#     with the assumption that a consistent time zone applies for the time 
#     period.

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

# Set region in
#
region_in <- 1

# Set up directories to read in and output data
#
era_interdir <- paste0("/projectnb/gislane/Zach/CAFE_ERA5_WBGT/InterDir/region", region_in, "/") # Directory where WBGT and UTCI are output
outdir <- paste0("/projectnb/gislane/Zach/CAFE_ERA5_WBGT/OutputData/region", region_in, "/")     # Output directory for county-aggregated WBGT/UTCI
points_dir <- "/projectnb/gislane/Zach/CAFE_ERA5_WBGT/RawData/ERA5_Fishnet/" # Input directory for extraction points (5A)

# Time zone specification. 
# The SpatRaster as downloaded from Copernicus will include hourly data based on
# the UTC time zone. When we calculate our daily summary statistics in the loop
# below, we want to make sure we are calculating statistics from midnight to
# midnight, using the local time zone. Specify the time zone relevant for your 
# data below. Note: If you are attempting to aggregate across multiple distinct
# time zones, additional processing is necessary. This capability is in development
# for an additional version of this process. To see a list of all time zones,
# run: OlsonNames()
#
tz_country <- "US/Eastern"

# Set name of geographic feature across which you want to summarize. This
# should be a column in the extraction_pts datasets
#
agg_geo <- "GEOID"

# Set year to process
#
years_to_agg <- c(2000:2004)

################### Run Processing to Daily Aggregate ##########################

# Define a function to calculate sums such that if all values are NA then it returns
# NA rather than 0.
#
sumfun <- function(x) { return(ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))) }

# Read in extraction points
#
extraction_pts <- vect(paste0(points_dir, "era5_county_usreg", region_in, "_extraction_pts.gpkg"))

# Time zone specification. 
# The SpatRaster as downloaded from Copernicus will include hourly data based on
# the UTC time zone. When we calculate our daily summary statistics in the loop
# below, we want to make sure we are calculating statistics from midnight to
# midnight, using the local time zone. To account for multiple time zones, we
# will identify all unique time zones in the spatial data using the lutz package.
# We will then run the daily aggregation as a loop of each unique time zone and
# stack the outputs (see line 144 to 320)
#
extraction_pts_sf <- st_as_sf(extraction_pts)
extraction_pts_sf$tz <- lutz::tz_lookup(extraction_pts_sf$geometry)
extraction_pts <- vect(extraction_pts_sf)
all_tz <- unique(extraction_pts$tz)

# Read in water (will be masked)
#
water_perm <- rast(paste0(era_interdir, "jrc_perm_era5.tif"))

# Set up path for all ERA5 wbgt and utci files
#
era_files <- list.files(path = era_interdir, pattern = ("utci.*.tif$|wbgt.*.tif$"))

for (year in c(years_to_agg)) {
  
  cat("Now processing ", year, "\n")
  
  # Subset to year from all files 
  #
  era_files_yr <- era_files[grepl(year, era_files) | grepl(year - 1, era_files) ]
  
  # Stack all of the daily files by year
  #
  era_files_yr <- paste0(era_interdir, "/", era_files_yr)
  
  # Read in UTCI and WBGT
  #
  wbgt <- rast(era_files_yr[grepl("wbgt", era_files_yr)])
  utci <- rast(era_files_yr[grepl("utci", era_files_yr)])
  
  # Mask out water tiles
  #
  wbgt <- mask(wbgt, water_perm == 1, inverse = T)
  utci <- mask(utci, water_perm == 1, inverse = T)
  
  # Filter through the unique time zones that are included
  #
  for (timezone in c(all_tz)) {
    
    ##################### Time Zone Conversion ###################################
    # Align times
    #
    terra::time(wbgt) <- terra::time(utci)
    
    # Reset times to align with specified time zone
    #
    terra::time(utci) <- with_tz(terra::time(utci) , tzone = timezone)
    terra::time(wbgt) <- with_tz(terra::time(wbgt) , tzone = timezone)
    
    # Subset stack to exclude the times that run past specified year due to 
    # time zone adjustment, and to exclude the previous year that was read in 
    # for time zone adjustment
    #
    utci <- subset(utci, time(utci) < date(paste0(year + 1, "-01-01")) &
                     time(utci) >= date(paste0(year, "-01-01")))
    wbgt <- subset(wbgt, time(wbgt) < date(paste0(year + 1, "-01-01")) &
                     time(wbgt) >= date(paste0(year, "-01-01")))
    
    # Get n layers
    #
    layer_n <- nlyr(utci)
    
    ##################### Daily Aggregation ######################################
    # Create a time sequence starting from January first of the year date
    #
    start_date <- as.POSIXct(paste0(year, "-01-01 00:00"), tz = timezone)
    time_seq <- seq(from = start_date, by = "hour", length.out = layer_n)
    
    # Convert our time sequence to a factor format. This will allow for use as 
    # a grouping variable in assessing daily level summary measures
    #
    daily_factor <- as.factor((substr(time_seq, 1, 10)))
    
    # Ensure time zone is assigned
    #
    time(utci) <- with_tz(terra::time(utci) , tzone = timezone)
    time(wbgt) <- with_tz(terra::time(wbgt) , tzone = timezone)
    
    # Set list of rasters. We will do the same processing of daily minimum, mean,
    # and maximum for the three variables we queried through ERA5, so we set
    # the lasters in a list and conduct the processing as below
    #
    list_rasters <- list(utci, wbgt)
    
    for (i in 1:length(list_rasters)) {
      
      # Tracker for viewing progress
      #
      cat("Now processing raster ", i, "\n")
      
      # Aggregate to daily mean temperature
      #
      daily_mean <- tapp(list_rasters[[i]], daily_factor, fun = mean)
      
      # Aggregate to daily maximum temperature
      #
      daily_max <- tapp(list_rasters[[i]], daily_factor, fun = max)
      
      # Aggregate to daily minimum temperature
      #
      daily_min <- tapp(list_rasters[[i]], daily_factor, fun = min)
      
      # Project points to WGS84 (coordinate system for ERA stack)
      #
      extraction_pts_cut <- extraction_pts[extraction_pts$tz == timezone, ]
      extraction_pts_cut <- project(extraction_pts_cut, crs(list_rasters[[i]]))
      
      # Extract daily summaries to point-based grid
      #
      mean_pts <- terra::extract(daily_mean, extraction_pts_cut)
      max_pts <- terra::extract(daily_max, extraction_pts_cut)
      min_pts <- terra::extract(daily_min, extraction_pts_cut)
      
      # Join results with extraction points (includes geographic identifiers)
      #
      mean_pts <- cbind(extraction_pts_cut, mean_pts)
      max_pts <- cbind(extraction_pts_cut, max_pts)
      min_pts <- cbind(extraction_pts_cut, min_pts)
      
      ####################### Link Missing Data ##################################
      # Note that here were are linking points that are missing to the nearest
      # value. This may not be advisable in all contexts! Exploring where missing
      # points are is important to determine whether leaving as missing or joining
      # to nearest available points is advisable
      #
      # Use package sf to join points on coast and not covered by ERA5 land to
      # nearby temperature measures
      #
      mean_pts_sf <- st_as_sf(mean_pts)
      max_pts_sf <- st_as_sf(max_pts)
      min_pts_sf <- st_as_sf(min_pts)
      
      # Identify a date in year to assess missingness for points outside
      # land extent
      #
      varname <- paste0("X", year, ".01.01")
      
      # Subset missing data
      #
      mean_pts_missing <- mean_pts_sf[is.na(mean_pts_sf[[varname]]), ]
      max_pts_missing <- max_pts_sf[is.na(max_pts_sf[[varname]]), ]
      min_pts_missing <- min_pts_sf[is.na(min_pts_sf[[varname]]), ]
      
      # Subset available data
      #
      mean_pts_avail <- mean_pts_sf[!is.na(mean_pts_sf[[varname]]), ]
      max_pts_avail <- max_pts_sf[!is.na(max_pts_sf[[varname]]), ]
      min_pts_avail <- min_pts_sf[!is.na(min_pts_sf[[varname]]), ]
      
      
      # Rejoin updated data with those already including temp
      #
      mean_pts <- rbind(mean_pts_avail, mean_pts_missing)
      max_pts <- rbind(max_pts_avail, max_pts_missing)
      min_pts <- rbind(min_pts_avail, min_pts_missing)
      
      # Convert the extracted data to a data frame
      #
      mean_pts_df <- as.data.frame(mean_pts)
      max_pts_df <- as.data.frame(max_pts)
      min_pts_df <- as.data.frame(min_pts)
      
      # Remove geometry column
      #
      mean_pts_df$geometry <- NULL
      max_pts_df$geometry <- NULL
      min_pts_df$geometry <- NULL
      
      # Extract columns relevant to ERA5 data
      #
      era5_cols <- names(mean_pts_df)[!names(mean_pts_df) %in% c("ID.1", names(extraction_pts))]
      
      # Set names for ERA5 variables based on input raster naming
      #
      mean_name <- paste0(substr(names(list_rasters[[i]])[1], 1, 3), "_mean")
      max_name <- paste0(substr(names(list_rasters[[i]])[1], 1, 3), "_max")
      min_name <- paste0(substr(names(list_rasters[[i]])[1], 1, 3), "_min")
      
      # Transpose the data frame to get time series format (long-form)
      #
      mean_long <- mean_pts_df %>%
        pivot_longer(cols = all_of(era5_cols), names_to = "date", values_to = mean_name) 
      
      max_long <- max_pts_df %>%
        pivot_longer(cols = all_of(era5_cols), names_to = "date", values_to = max_name) %>%
        select(UniqueID, date, !!sym(max_name))
      
      min_long <- min_pts_df %>%
        pivot_longer(cols = all_of(era5_cols), names_to = "date", values_to = min_name) %>%
        select(UniqueID, date, !!sym(min_name))
      
      # Combine data into single dataframe
      #
      era5_long <- left_join(mean_long, max_long, by = c("UniqueID", "date")) %>%
        left_join(., min_long, by = c("UniqueID", "date"))
      
      # Save data as date format
      #
      era5_long$date <- as.Date(substr(era5_long$date, 2, 11), format = "%Y.%m.%d")
      
      # Join together all measures
      #
      if (i == 1) {
        era5_full <- era5_long
      } else if ( i != 1 ) {
        era5_long <- era5_long[c("UniqueID", "date", mean_name, max_name, min_name)]
        era5_full <- left_join(era5_full, era5_long, by = c("UniqueID", "date"))
      }
      
    }
    
    # Stack to create full dataset with all time zones
    #
    if (timezone == all_tz[1]) {
      era5_complete <- era5_full
    } else {
      era5_complete <- rbind(era5_complete, era5_full)
    }
  }
  
  
  
  ############################## County Aggregation ############################
  # In this step, we will use the extraction points to extract the ERA5 grid cell
  # underlying each portion of a ward across the entire raster stack of values.
  # We again will follow the same process for each individual variable, and use 
  # a loop to conduct the processing
  #
  varnames <- names(era5_complete)[!names(era5_complete) %in% c(names(extraction_pts), "date", "ID.1")]
  
  for (i in 1:length(varnames)) {
    
    cat("Now processing ", varnames[i], "\n")
    
    # Before we calculate the final weighted average of the ERA5 measure, we need to check for missing data.
    # If a value is NA on one of the polygons, then it will give an underestimate of the
    # temperature since the weights will no longer add to 1. Example: there are two
    # polygons, each with 50% area. If Tmax is 30 C in one and NA in the other, then
    # the area weighted average (which removes NA values) would give: (30 * 0.5) + (NA * 0.5) = 15 C.
    # Therefore, we need to re-weight the weights based on the availability of data.
    #
    eqn <- as.formula(paste0("SpatWt ~ ", agg_geo, " + date")) 
    avail <- summaryBy(eqn,
                       data = era5_complete[which( !(is.na(era5_complete[[varnames[i]]])) ),],
                       FUN = sumfun)
    
    # Merge this value back into the longform ERA5 data
    #
    era5_complete <- merge(era5_complete, avail, by = c(agg_geo, "date"), all.x = TRUE)
    
    # Re-weight the area weight by dividing by total available weight
    #
    era5_complete$SpatWt <- era5_complete$SpatWt / era5_complete$SpatWt.sumfun
    
    # QC: check that the weights of *available data* all add to 1
    #
    eqn <- as.formula(paste0("SpatWt ~ ", agg_geo, " + date"))
    check <- summaryBy(eqn,
                       data = era5_complete[which( !(is.na(era5_complete[[varnames[i]]])) ),],
                       FUN = sumfun)
    
    if (length(which(round(check$SpatWt.sumfun, 4) != 1)) > 0) {
      cat("ERROR: weights do not sum to 1", "\n"); break 
    } else {
      cat(":) weights sum to 1", "\n")
      era5_complete$SpatWt.sumfun <- NULL
    }
    
    # Multiply the variable of interest (here "newvarname") by the weighting value and then
    # sum up the resultant values within admin boundaries. This is an area-weighted average.
    #
    tempvar <- paste0(varnames[i], "_Wt")
    era5_complete[[tempvar]] <- era5_complete[[varnames[i]]] * era5_complete[["SpatWt"]]
    
    eqn <- as.formula(paste0(tempvar, " ~ ", agg_geo, " + date")) 
    
    final <- summaryBy(eqn, data = era5_complete, FUN = sumfun)
    
    # Automated QC to confirm that the dimensions are correct
    #
    if ( length(unique(as.data.frame(extraction_pts)[[agg_geo]]))  * length(unique(era5_complete$date)) != dim(final)[1]) {
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
  missing_utci_max <- which(is.na(finaloutput$utci_mean))
  missing_utci_min <- which(is.na(finaloutput$utci_min))
  missing_utci_mea <- which(is.na(finaloutput$utci_max))
  missing_wbgt_max <- which(is.na(finaloutput$wbgt_mean))
  missing_wbgt_min <- which(is.na(finaloutput$wbgt_min))
  missing_wbgt_mea <- which(is.na(finaloutput$wbgt_max))
  
  if (length(missing_utci_max) > 0 | length(missing_wbgt_max) > 0) {
    cat("WARNING: Note the number of missing ward-days by variable: \n")
    cat("UTCI Max:", length(missing_utci_max), "\n")
    cat("WBGT Max:", length(missing_wbgt_max), "\n")
    
    cat("The first few lines of missing UTCI (if any) are printed below: \n")
    print(head(finaloutput[missing_utci_max,]))
    
    cat("The first few lines of missing WBGT (if any) are printed below: \n")
    print(head(finaloutput[missing_wbgt_max,]))
    
  } else { cat(":) No missing temperature values! \n") }
  
  # Automated QC: impossible temperature values
  #
  num_temp_errors_utci <- length(which(finaloutput$utci_max < finaloutput$utci_mean |
                                         finaloutput$utci_max < finaloutput$utci_min |
                                         finaloutput$utci_min > finaloutput$utci_mean |
                                         finaloutput$utci_min > finaloutput$utci_max))
  
  num_temp_errors_wbgt <- length(which(finaloutput$wbgt_max < finaloutput$wbgt_mean |
                                         finaloutput$wbgt_max < finaloutput$wbgt_min |
                                         finaloutput$wbgt_min > finaloutput$wbgt_mean |
                                         finaloutput$wbgt_min > finaloutput$wbgt_max))
  
  if (num_temp_errors_utci > 0 | num_temp_errors_wbgt > 0 ) { 
    
    print("ERROR: impossible temperature values. Applicable rows printed below:")
    print(finaloutput[which(finaloutput$utci_max < finaloutput$utci_mean |
                              finaloutput$utci_max < finaloutput$utci_min |
                              finaloutput$utci_min > finaloutput$utci_mean |
                              finaloutput$utci_min > finaloutput$utci_max),])
    print(finaloutput[which(finaloutput$wbgt_max < finaloutput$wbgt_mean |
                              finaloutput$wbgt_max < finaloutput$wbgt_min |
                              finaloutput$wbgt_min > finaloutput$wbgt_mean |
                              finaloutput$wbgt_min > finaloutput$wbgt_max),])
    
  } else { print(":) all temperature values are of correct *relative* magnitude") }
  
  # Output results by year to output directory
  #
  saveRDS(finaloutput, paste0(outdir, "/", "county_agg_era5_", year, "_wbgt_utci_region", region_in, ".rds"))
  
}

