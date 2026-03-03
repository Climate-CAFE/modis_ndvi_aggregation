# Created by: Zach Popp
# Adapted by: Allison James
# Date Created: 02/19/2026
# Version Number: v1
# Date Modified: 
# Modifications:
#   Code modified from: https://github.com/Climate-CAFE/ERA5-daily-heat-aggregation
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
#     first in a two-step raster processing process (following the steps to 
#     create hourly WBGT and UTCI in scripts 1-4). In this script we build
#     a polygon grid that intersects our polygons of interest against the grid
#     of MODIS WBGT and UTCI inputs. This will be used to calculate area-weighted
#     estimates of heat metrics at the polygon level

# Load required packages
#
library("terra")  # For raster data
library("sf")     # For vector data
library("plyr")
library("doBy")
library("tidyverse")
library("tidycensus")
library("lwgeom")
library("tigris")

sf_use_s2(FALSE)  
# S2 is for computing distances, areas, etc. on a SPHERE (using
# geographic coordinates, i.e., lat/lon in decimal-degrees); no need for this
# extra computational processing time if using PROJECTED coordinates,
# since these are already mapped to a flat surface. Here, MODIS data
# is indeed in geographic coordinates, but the scale of areas we are 
# interested in is very small, and hence the error introduced by 
# ignoring the Earth's curvature over these tiny areas is negligible and
# a reasonable trade off given the dramatic reduction in processing time. Moreover,
# the areas we calculate are not an integral part of the process
# and any error in that step would not materially impact the final output

# Check package version numbers
#
if (packageVersion("terra") < "1.5.34"   | packageVersion("sf") < "1.0.7" | 
    packageVersion("plyr")  < "1.8.7"    | packageVersion("lwgeom") < "0.2.8" | 
    packageVersion("doBy")  < "4.6.19"   | packageVersion("tidyverse") < "1.3.1") {
  cat("WARNING: packages are outdated and may result in errors.") }

########################## User-Defined Parameters #############################


# Set up directories to read in and output data
# 

dir <- "/Users/alliej/Documents/cafe/modis_ndvi_data/"
modis_dir <- paste0(dir, "Suffolk/summer_mean/")
points_dir <- paste0(dir, "MODIS_Fishnet/")
  
  
# LOAD Shapefile. This approach involves a US application for Suffolk County, 
# Massachusetts downloaded using TIGRIS. If you are conducting a global analysis or have 
# an existing shapefile, the below can be replaced to just set the input_shape
# input (this script does not use the shapefile as before, it just uses US 
# state names from that shapefle to download US counties)
#     input_shape <- st_read("shapefile_path")
#
#shapefile <- tigris::counties(year = 2020)

# Subset to the actual region of interest
#
region_in <- "25025"
#shapefile <- shapefile[shapefile$GEOID == region_in, ]

# Query tracts in region - we first use the above to set the counties we need,
# then download the actual shape here. This approach could be easily modified
# to query counties in a state, for example.
#
state <- "MA" #Feel free to delete or change
county_name <- "Suffolk" #Feel free to delete or change
input_shape <- tigris::tracts(state = state, county = county_name, year = 2020)

# Remove water - since we have masked out water in the MODIS data already,
# we don't want water to be included in NDVI aggregations later. For TIGER/Line
# shapefiles, we can just remove polygons that have a land area of 0.
#
input_shape <- input_shape %>% filter(ALAND != 0)

####################### Develop Grid ###########################################

# Read in MODIS raster data
#
modis_files <- list.files(modis_dir, pattern=paste0('mean.*.tif'), full.names = F)

# Stack all of the files
#
modis_files <- paste0(modis_dir, "/", modis_files)

# We just need one layer, as the grid is the same across the board
#
modis_stack <- rast(modis_files[1])

# %%%%%%%%%%%%%%%%%%%% CREATE A FISHNET GRID OF THE RASTER EXTENT %%%%%%%%%%%% #
#
# Here, we are making a shapefile that is a fishnet grid of the raster extent.
# It will essentially be a polygon of lines surrounding each MODIS cell.
#
# Reference/credit: https://gis.stackexchange.com/a/243585
#
modis_raster <- modis_stack[[1]]
modis_extent <- ext(modis_raster)
xmin <- modis_extent[1]
xmax <- modis_extent[2]
ymin <- modis_extent[3]
ymax <- modis_extent[4]

modis_matrix <- matrix(c(xmin, ymax,
                        xmax, ymax,
                        xmax, ymin,
                        xmin, ymin,
                        xmin, ymax), byrow = TRUE, ncol = 2) %>%
  list() %>% 
  st_polygon() %>% 
  st_sfc(., crs = st_crs(modis_raster))

# Create fishnet of the MODIS matrix. This takes some time.
#
modis_rows <- dim(modis_raster)[1]
modis_cols <- dim(modis_raster)[2]
modis_fishnet <- st_make_grid(modis_matrix, n = c(modis_cols, modis_rows), 
                             crs = st_crs(modis_raster), what = 'polygons') %>%
  st_sf('geometry' = ., data.frame('ID' = 1:length(.)))


# Automated QC check -- confirm same coordinate reference system (CRS) between
#                       the fishnet and MODIS raster
if ( !(all.equal(st_crs(modis_raster), st_crs(modis_fishnet))) ) {
  cat("ERROR: CRS's do not match \n") } else { cat(":) CRS's match \n") }

############# Form Fishnet - Shape Union #######################################

# Rename geom column to geometry if necessary, to align with fishnet
#
if (!("geometry" %in% names(input_shape))) {
  input_shape <- input_shape %>% rename(geometry = geom)
}

# Match the CRS of the wards shapefile to the fishnet and era data and confirm match
#
input_shape <- st_transform(input_shape, crs = st_crs(modis_fishnet))

# Run check to ensure CRS are equal
#
if (!isTRUE(all.equal(st_crs(input_shape), st_crs(modis_fishnet)))) {
  cat("ERROR: CRS's don't match \n")  } else { cat(":) CRS's match \n") }

# %%%%%%%%%%%%%%%%%% CREATE UNION BETWEEN FISHNET AND WARDS %%%%%%%%% #
# Reference/credit: https://stackoverflow.com/a/68713743
#
my_union <- function(a,b) {
  st_agr(a) = "constant"
  st_agr(b) = "constant"
  op1 <- st_difference(a, st_union(b))
  op2 <- st_difference(b, st_union(a))
  op3 <- st_intersection(b, a)
  union <- rbind.fill(op1, op2, op3)
  return(st_as_sf(union))
}

# Ensure geometries are valid
#
input_shape <- st_make_valid(input_shape) 
modis_fishnet <- st_make_valid(modis_fishnet)

# Create the union between the fishnet and blocks layers
#
fishnetward <- my_union(modis_fishnet, input_shape)
fishnetward$UniqueID <- 1:dim(fishnetward)[1]

# Automated QC -- Check to see if the union has introduced any geometry errors
#                 and fix as appropriate
#
check <- try(st_make_valid(fishnetward), silent = TRUE)

if (class(check)[1] == "try-error") {
  
  cat("There is an issue with the sf object \n")
  cat("..... Attempting fix \n")
  
  geo_types <- unique(as.character(st_geometry_type(fishnetward, by_geometry = TRUE)))
  
  cat("..... Geometry types in sf object 'fishnetward':", geo_types, "\n")
  
  for (j in 1:length(geo_types)) {
    fishnetward_subset <- fishnetward[which(st_geometry_type(fishnetward, by_geometry = TRUE) == geo_types[j]),]
    if (j == 1) { updated_fishnetward <- fishnetward_subset; next }
    updated_fishnetward <- rbind(updated_fishnetward, fishnetward_subset)
    
  }
  
  check2 <- try(st_make_valid(updated_fishnetward), silent = TRUE)
  
  if (class(check2)[1] == "try-error") {
    cat("..... ERROR NOT RESOLVED \n") } else {
      cat("..... :) issue has been fixed! \n")
      
      updated_fishnetward <- updated_fishnetward[order(updated_fishnetward$UniqueID),]
      if ( !(all.equal(updated_fishnetward$UniqueID, fishnetward$UniqueID)) ) {
        cat("ERROR: Unique IDs do not match \n") } else {
          cat(":) unique ID's match. Reassigning 'updated_fishnetward' to 'fishnetward' \n")
          fishnetward <- updated_fishnetward    
        }
    }
}

# Automated QC check -- ensure that there is a variable identifying the geographies
#                       that will be used for aggregation (here tracts). 
#                       Note that these variable names may change
#                       depending on the version and country used. 
# Specify name of the geographic identifier in your administrative boundary data
#
geo_name <- "GEOID"

# Extract name from dataset
#
geo_id_var <- names(fishnetward)[grep(paste0("^", geo_name), names(fishnetward), ignore.case = TRUE)]

# Identify the polygons of the fishnet that do not intersect with the ward
# data; drop them.
#
before_dim <- dim(fishnetward)[1]
fishnetward <- fishnetward[which( !(is.na(fishnetward[[geo_id_var]])) ),]
after_dim <- dim(fishnetward)[1]

cat("Dropped", before_dim - after_dim, "polygons that do not intersect with census data \n")

# Some polygons formed in the union are incredibly small -- this adds unnecessary
# computation time without materially reducing error. Drop the small polygons.
# NOTE: Typically, when calculating areas of polygons, you would want to convert to
#       a projected CRS appropriate for your study domain. For the purpose of identifying
#       negligibly small areas to drop here, the error introduced by using geographic
#       coordinates for calculating area at this scale is negligible.
#
fishnetward$Area_m2 <- as.numeric(st_area(fishnetward))

fishnetward <- fishnetward[which(fishnetward$Area_m2 > 10),]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERT POLYGON TO POINTS %%%%%%%%%%%%%%%%%%%%% #
#
# The final step is to create the extraction points. This is a point shapefile
# that will enable us to extract MODIS data from an entire stack of rasters rather
# than individually processing zonal statistics on each raster layer. 
#
# NOTE: This step throws a warning message related to using geographic coordinates 
#       rather than a projected CRS. This step is only placing a point inside the 
#       polygon to identify which MODIS grid cell we need to extract from; as all
#       of the input data are on the same CRS and the spatial scale of the polygons
#       is extremely small, this does not introduce substantive error.
#
extraction_pts <- st_point_on_surface(fishnetward)

# Get the total area by ward to calculate the spatial weight value (typically 1.0)
#
eqn <- as.formula(paste0("Area_m2 ~ ", geo_id_var))

# Define a function to calculate sums such that if all values are NA then it returns
# NA rather than 0.
#
sumfun <- function(x) { return(ifelse(all(is.na(x)), NA, sum(x, na.rm = TRUE))) }

# Calculate the sum area by ward
#
ptstotal <- summaryBy(eqn, data = as.data.frame(extraction_pts), FUN = sumfun)

# Merge area and calculate spatial weight of points
#
extraction_pts <- merge(extraction_pts, ptstotal, by = geo_id_var, all.x = TRUE)
extraction_pts$SpatWt <- extraction_pts$Area_m2 / extraction_pts$Area_m2.sumfun

# Save extraction points
#
st_write(extraction_pts, paste0(points_dir, "MODIS_county_usreg", region_in, "_extraction_pts.gpkg"),
         append = FALSE)
