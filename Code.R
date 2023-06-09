# Daniel Lauer
# 26 May 2021 - 9 June 2023

# Paper title: Africa's ecosystems exhibit a tradeoff between resistance and stability following disturbances

###########################################################################################################################
###########################################################################################################################

# PART I: SELECTION OF SITES IN AFRICAN PROTECTED AREAS

# Set the working directory to the folder containing the relevant data sources:

setwd('./Data')

# Obtain and plot polygons of all African countries, taken from a global map provided by Natural Earth:

library(raster) # For using the "shapefile" function, and for general spatial analyses throughout.
Africa <- shapefile('Shapefile_Africa/Africa.shp') # Shapefile of whole planet.
Africa <- Africa[Africa@data$CONTINENT == 'Africa', ] # Subset to just Africa.
plot(Africa, col = 'cornsilk1') # Plot to ensure it looks like Africa.

# Obtain, subset, and process polygons of global WDPA protected areas (PAs), which are the study units for this chapter
# (elements of this code block were run once and commented out, with their outputs saved, to save time and space):

  # # Obtain the PAs from the WDPA database as a list of the three shapefiles provided by the WDPA:
  # 
  # Africa_PAs <- list() # List to hold the PA polygons that occur within "Africa".
  # for (File in 1:3) { # For each of the three shapefiles provided by the WDPA...
  #   Africa_PAs[[File]] <- shapefile(paste('./WDPA_May2021_Public_shp_', as.character(File-1),
  #     '/WDPA_May2021_Public_shp-polygons.shp', sep = '')) # Read it in.
  #   print(paste0('Shapefile ', File, ' Read In')) }; rm(File) # Keep a record of this loop's progress.
  # 
  # # Subset each WDPA shapefile down to just the PAs that are terrestrial and occurring within the bounds of "Africa":
  # 
  # Africa_PAs <- lapply(Africa_PAs, function(Item) { # For each shapefile in "Africa_PAs"...
  #   Item <- Item[Item@data$MARINE != '2', ] # Subset to just terrestrial PAs.
  #   Item <- crop(Item, Africa) }) # Crop to just PAs in "Africa".
  # 
  # # Compress the WDPA shapefiles to a single layer and subset the layer to just the PAs established before 2000:
  # 
  # Africa_PAs <- bind(Africa_PAs[[1]], Africa_PAs[[2]], Africa_PAs[[3]]) # Compress "Africa_PAs" into a single layer.
  # library(rgdal) # For using the "writeOGR" function.
  # writeOGR(Africa_PAs, 'Shapefile_AfricaPAs', layer = 'AfricaPAs', driver = 'ESRI Shapefile') # Save "Africa_PAs".
  Africa_PAs <- shapefile('Shapefile_AfricaPAs/AfricaPAs.shp') # Read in "Africa_PAs".
  Africa_PAs <- Africa_PAs[Africa_PAs@data$STATUS_YR < 2000, ] # Subset based on year.
  # 
  # # Obtain the centroid of each remaining PA as a point:
  # 
  library(rgeos) # For using the "gCentroid" and "gIntersection" functions.
  Africa_PAs_centroids <- gCentroid(Africa_PAs, byid = TRUE) # Calculate the centroids.
  # 
  # # Loop through the polygons in "Africa_PAs". Per polygon, identify the raster pixel of a given resolution (the
  # # largest resolution to be considered in this chapter) that surrounds the polygon's centroid. Convert the pixel into
  # # a polygon, and determine if the original polygon completely contains the pixel. Then, subset the PAs in question
  # # down to just those that achieve this complete containment:
  # 
  # Africa_PAs_rasters <- list(raster(list.files()[grep('1km_uint', list.files())]), raster(list.files()[grep('5km_uint',
  #   list.files())])) # Examples of 1x1-km and 5x5-km global rasters, taken from EarthEnv.
  # Africa_PAs_rasters <- lapply(Africa_PAs_rasters, function(Ras) { # For each raster at each resolution...
  #   raster(crs = crs(Ras), ext = extent(Africa), resolution = res(Ras)) }) # Blank it out and crop it to just "Africa".
  # 
  # Africa_PAs_contains <- c() # Vector to hold if each polygon achieves complete containment.
  # for (PA in 1:nrow(Africa_PAs)) { # For each PA polygon in "Africa_PAs"...
  #   Tmp <- rasterize(Africa_PAs_centroids[PA, ], Africa_PAs_rasters[[length(Africa_PAs_rasters)]], field = 1) # ID pixel.
  #   Tmp <- rasterToPolygons(Tmp) # Convert pixel to polygon.
  #   Tmp <- spTransform(Tmp, crs(Africa_PAs)) # Re-project pixel polygon into the same CRS as "Africa_PAs".
  #   Tmp2 <- gIntersection(Tmp, Africa_PAs[PA, ]) # Get intersection between pixel polygon and PA polygon.
  #   Africa_PAs_contains[PA] <- ifelse(!is.null(Tmp2), area(Tmp2)/area(Tmp) >= 0.99, FALSE) # Complete containment?
  #   print(paste0('PA ', PA, ' complete')) }; rm(PA, Tmp, Tmp2) # Keep a record of this loop's progress.
  # save(Africa_PAs_contains, file = 'Africa_PAs_contains.rda') # Save "Africa_PAs_contains".
  load('Africa_PAs_contains.rda') # Read in "Africa_PAs_contains".
  Africa_PAs <- Africa_PAs[Africa_PAs_contains, ]; Africa_PAs_centroids <- Africa_PAs_centroids[Africa_PAs_contains, ]
  # 
  # # Identify the raster pixel of all resolutions of interest that surrounds each remaining PA's centroid. Some pixels,
  # # especially at the coarsest resolution, may contain more than one centroid, due to the overlap of certain PAs.
  # # In those cases, keep only the PA with the largest area, and remove the others:
  # 
  # Africa_PAs_rasters <- lapply(Africa_PAs_rasters, function(Ras) { # For each raster at each resolution...
  #   rasterize(Africa_PAs_centroids, Ras, fun = 'count') }) # Identify pixels, and count number of centroids in each.
  # 
  # Africa_PAs_overlap_pix <- which(values(Africa_PAs_rasters[[length(Africa_PAs_rasters)]]) > 1) # Pixels with overlap.
  # Africa_PAs_overlaps <- list() # List to hold which overlapping PAs should be removed.
  # Counter <- 1 # Counter for the loop below.
  # 
  # for (Pixel in Africa_PAs_overlap_pix) { # For each pixel with PA overlap...
  #   Tmp <- Africa_PAs_rasters[[length(Africa_PAs_rasters)]] # Copy the raster that contains the pixel.
  #   values(Tmp) <- NA; values(Tmp)[Pixel] <- 1 # Make all raster values NA, except for that of the pixel.
  #   Africa_PAs_overlaps[[Counter]] <- which(!is.na(extract(Tmp, Africa_PAs_centroids))) # Which PAs in that pixel.
  #   Africa_PAs_overlaps[[Counter]] <- Africa_PAs_overlaps[[Counter]][-which.max(Africa_PAs@data$GIS_AREA[
  #     Africa_PAs_overlaps[[Counter]]])] # Keep a record of the smaller PAs that will be removed.
  #   Counter <- Counter + 1 }; rm(Counter, Pixel, Tmp) # Update the counter.
  # 
  # Africa_PAs_overlaps <- unlist(Africa_PAs_overlaps) # Make a vector of the PAs to be removed.
  # save(Africa_PAs_overlaps, file = 'Africa_PAs_overlaps.rda') # Save "Africa_PAs_overlaps".
  load('Africa_PAs_overlaps.rda') # Read in "Africa_PAs_overlaps".
  Africa_PAs <- Africa_PAs[-Africa_PAs_overlaps, ]; Africa_PAs_centroids <- Africa_PAs_centroids[-Africa_PAs_overlaps, ]
  # 
  # Africa_PAs_rasters <- lapply(Africa_PAs_rasters, function(Ras) { # For each raster at each resolution...
  #   rasterize(Africa_PAs_centroids, Ras, field = Africa_PAs_centroids@coords[,2]) }) # Put centroid latitudes in pixels.
  # save(Africa_PAs_rasters, file = 'Africa_PAs_rasters.rda') # Save "Africa_PAs_rasters".
  load('Africa_PAs_rasters.rda') # Read in "Africa_PAs_rasters".
  # 
  # # Plot the centroids of the remaining PAs on top of "Africa". Add a scale bar to the plot:
  # 
  library(maps) # For using the "map.scale" function.
  plot(Africa_PAs_centroids, col = 'cornflowerblue', pch = 19, add = TRUE) # Plot centroids on "Africa".
  map.scale(ratio = FALSE, cex = 1.5) # Scale bar.
  
###########################################################################################################################
###########################################################################################################################

# PART II: COLLECTION OF DATA ADDRESSING THE ATTRIBUTES OF THE PROTECTED AREA SITES
  
# Per spatial resolution of interest, collect data on each PA site that is provided by or derivable from WDPA data
# (elements of this code block are commented out, as they are saved and read in anyways below):

Sites_res <- c(1,5) # The spatial resolutions of interest in XxX-km. For reference below.
# library(geosphere) # For using the "areaPolygon" and "perimeter" functions.
# 
# Sites_desigs <- Africa_PAs@data$DESIG_ENG # Initial designations for each PA, of which there are many.
# Sites_designator <- function(Desigs, Keywords, Name) { # Create a function that consolidates PA designations.
#   Desigs[grep(Keywords, Desigs, ignore.case = TRUE)] <- Name # Consolidate PAs based on designation keywords.
#   return(Desigs) } # Return the result after consolidation.
# Sites_desigs_pairs <- list(c('Hunt|Game', 'Hunting Area'), c('Classified', 'Classified Forest'), c(
#   'Bird|Fauna|Wildlife', 'Wildlife Reserve'), c('World|UNESCO|International', 'International Site'),
#   c('Absolute|Strict|Integral|Priority|Wilderness', 'Strict Reserve'), c('Nature Reserve|Natural Reserve',
#   'Nature Reserve')) # Keywords paired with full names of new PA designations, for consolidation.
# for (Pair in Sites_desigs_pairs) { # Per keyword-designation name pair...
#   Sites_desigs <- Sites_designator(Sites_desigs, Pair[1], Pair[2]) }; rm(Pair) # Perform the consolidation using pair.
# Sites_desigs[!(Sites_desigs %in% c(sapply(Sites_desigs_pairs, '[[', 2), 'Forest Reserve', 'National Park'))] <- 'Other'
# Sites_desigs <- as.factor(Sites_desigs) # Convert to factor.
# 
# Sites_IUCN <- Africa_PAs@data$IUCN_CAT # Initial IUCN categories for each PA.
# Sites_IUCN[grep('Not', Sites_IUCN)] <- 'Unrecorded' # Consolidate IUCN categories for which one was not recorded.
# Sites_IUCN <- as.factor(Sites_IUCN) # Convert to factor.
# 
# Sites <- replicate(length(Sites_res), data.frame(Site_Lat = Africa_PAs_centroids@coords[,2], Site_Long =
#   Africa_PAs_centroids@coords[,1], Site_Country = Africa_PAs@data$ISO3, Site_ID = Africa_PAs@data$WDPA_PID,
#   PA_Area = areaPolygon(Africa_PAs)/1000, PA_Desig = Sites_desigs, PA_IUCN = Sites_IUCN), simplify = FALSE) # Df per res.
# Sites <- lapply(Sites, transform, # For each dataframe, one per spatial resolution...
#   PA_Shape = sqrt(PA_Area / (pi * ((2 * (perimeter(Africa_PAs)/1000)) / pi)^2))) # Input PA shape index in each.

# Per spatial resolution, collect data on each PA site that involves intensive spatial data integration (this code block
# was run once per data variable collected and commented out, with its outputs saved each time, to save time and space):

  # # Read in the prior version of "Sites", which will be referenced and added to below:
  # 
  # load('Sites_V2_OneSigma+Full.rda')
  # 
  # # Collect data for all variables taken from single raster layers (uncomment and process each variable one by one):
  # 
  # # Var_name <- 'Plant_LDMC' # Variable to address each site's plants' leaf dry matter content (LDMC).
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LDMC_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_LDMCsd' # Variable to address each site's plants' standard deviation of LDMC.
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LDMC_sd_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_LNC' # Variable to address each site's plants' leaf nitrogen content (LNC).
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LNC_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_LNCsd' # Variable to address each site's plants' standard deviation of LNC.
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LNC_sd_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_LPC' # Variable to address each site's plants' leaf phosphorus content (LPC).
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LPC_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_LPCsd' # Variable to address each site's plants' standard deviation of LPC.
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/LPC_sd_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_SLA' # Variable to address each site's plants' specific leaf area (SLA).
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/SLA_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Plant_SLAsd' # Variable to address each site's plants' standard deviation of SLA.
  # # Var_data <- raster('Global_trait_maps_Moreno_Martinez_2018_Version2_1km_resolution/SLA_sd_1km_v1.tif') # Raster.
  # # Var_res <- 1 # Resolution for the variable. It is slightly different than 1-km "Africa_PAs_rasters", so re-project.
  # 
  # # Var_name <- 'Env_HetHab' # Variable to address each site's habitat heterogeneity (coefficient of variation in EVI).
  # # Var_data <- raster('cv_01_05_1km_uint16.tif') # Raster for the variable.
  # # Var_res <- 1 # Resolution for the variable.
  # 
  # # Var_name <- 'Env_HetTopo' # Variable to address each site's topographic heterogeneity (TRI).
  # # Var_data <- raster('tri_1KMmn_GMTEDmd.tif') # Raster for the variable.
  # # Var_res <- 1 # Resolution for the variable.
  # 
  # # Var_name <- 'Hum_TravelTime' # Variable to address each site's travel time to the nearest major city.
  # # Var_data <- raster('access_50k/acc_50k.tif') # Raster for the variable.
  # # Var_res <- 1 # Resolution for the variable.
  # 
  # Var_name <- 'Hum_HII' # Variable to address each site's Human Influence Index, or HII.
  # Var_data <- raster('hii-africa-geo-grid/hii_africa_grid/hii_africa/w001001.adf') # Raster for the variable.
  # Var_res <- 1 # Resolution for the variable.
  # 
  # if (!(identical(crs(Var_data), crs(Africa_PAs_rasters[[which(Sites_res == Var_res)]])))) { # If raster in wrong CRS...
  #   Var_data <- projectRaster(Var_data, Africa_PAs_rasters[[which(Sites_res == Var_res)]]) } else { # Re-project.
  #   Var_data <- crop(Var_data, Africa_PAs_rasters[[which(Sites_res == Var_res)]]) } # Or else, just crop to Africa.
  # Var_data <- unlist(list(Var_data, sapply(Sites_res[Sites_res > Var_res], function(Res) { aggregate(Var_data, Res,
  #   expand = FALSE) }))) # Add rasters to the variable representing coarser spatial resolution(s) of interest.
  # for (Res in seq(length(Sites_res))) { # For each resolution of interest...
  #   Sites[[Res]][Var_name] <- extract(Var_data[[Res]], Africa_PAs_centroids) }; rm(Res) # Add variable to res's df.
  # 
  # # Collect data for the variables representing ecosystem stability, which are taken from a stack of raster layers:
  # 
  #   # Read in the necessary outside function for this procedure, and define the names of the variables:
  # 
  #   source('../PhenoDeriv.R') # From the "greenbrown" R package, to calculate growing season timings.
  #   Var_name <- c('Eco_Resistance', 'Eco_RecoveryTime', 'Eco_RecoveryRate', 'Eco_Invariability') # Stability variables.
  # 
  #   # Process a series of URLs, where each URL points to a raster of EVI values, at 1-km resolution surrounding Africa,
  #   # for a specific 16-day period in a specific year throughout the time series of interest. The URLs were acquired
  #   # using Javascript code on the HTML file where the URLs are located:
  # 
  #   Var_URLs <- read.table('EVI_URLs.txt', sep = ',', header = FALSE) # Contains list of comma-separated URLs.
  #   Var_URLs <- as.character(unname(unlist(Var_URLs[1,]))) # Convert list from dataframe to vector.
  #   Var_URLs <- paste('/vsicurl/', Var_URLs, sep = '') # Give URLs GDAL capability.
  #   Var_URLs <- rev(Var_URLs) # Reverse the order of the URLs, so they go from oldest year to most recent.
  # 
  #   # Create a vector of the year that each URL refers to:
  # 
  #   Var_years <- regmatches(Var_URLs, regexec('doy\\s*(.*?)\\s*_', Var_URLs)) # Extract portion of URLs referencing year.
  #   Var_years <- sapply(Var_years, '[[', 2) # Perform further extraction.
  #   Var_years <- substr(Var_years, 1, 4) # Obtain years, such that every URL is matched to a year.
  #   Var_years <- Var_years[-length(Var_years)] # Remove the last year "2020", as it only has a single URL.
  # 
  #   # Create a vector of the 16-day period index within each year that each URL refers to:
  # 
  #   Var_periods <- max(sapply(unique(Var_years), function(Yr) { sum(Var_years == Yr) })) # Max periods for a given year.
  #   Var_periods <- unname(unlist(sapply(unique(Var_years), function(Yr) { # Per year...
  #     return((Var_periods - sum(Var_years == Yr) + 1):Var_periods) }))) # Return the index of each URL.
  # 
  #   # Create a raster stack of all rasters, represented by all URLs processed above:
  # 
  #   Var_data <- stack(Var_URLs[1:2]) # Create raster stack of first two rasters from first two URLs.
  #   for (Ras in 3:length(Var_URLs)) { # For the remaining URLs...
  #     Var_data <- addLayer(Var_data, Var_URLs[Ras]) # Add their rasters to the stack.
  #     print(paste('Raster', as.character(Ras), 'Stacked')) }; rm(Ras) # Keep a record of this loop's progress.
  #   save(Var_data, file = 'EVI_stack.rda') # Save "Var_data", which as of now contains a single raster stack.
  #   load('EVI_stack.rda') # Read in "Var_data".
  # 
  #   # Expand "Var_data" to include a raster stack not only for the native 1-km resolution, but for the other(s):
  # 
  #   Var_res <- Counter <- 1 # Resolution of the data provided by NASA in "Var_data", as well as a loop counter.
  #   Var_data <- list(Var_data) # Make the current version of "Var_data" the first element of a list.
  # 
  #   for (Res in Sites_res[Sites_res > Var_res]) { # For the remaining resolution(s) of interest...
  #     for (Layer in 1:length(Var_URLs)) { # For each raster layer in the current version of "Var_data"...
  #       Tmp <- aggregate(subset(Var_data[[1]], Layer), Res) # Reduce the layer's resolution to that of "Res".
  #       if (Layer == 1) { # If the first raster layer is being addressed...
  #         Var_data[[Counter+1]] <- Tmp } else { # Add reduced layer to a new list element of "Var_data". Otherwise...
  #         Var_data[[Counter+1]] <- addLayer(Var_data[[Counter+1]], Tmp) } # Stack the reduced layer on those previous.
  #       print(paste('Layer', as.character(Layer), 'Processed for', as.character(Res), 'km')) } # Progress record.
  #     Counter <- Counter + 1 }; rm(Res, Layer, Tmp, Counter) # Update the counter.
  # 
  #   save(Var_data, file = 'EVI_stacks.rda') # Save "Var_data", now that it has a stack for all resolutions.
  #   load('EVI_stacks.rda') # Read in "Var_data".
  # 
  #   # Per spatial resolution, create a dataframe in which each row refers to a site and each column to a 16-day period
  #   # in which EVI was measured. Fill the dataframe with the respective extracted values from "Var_data" for each site:
  # 
  #   Var_df <- lapply(1:length(Sites_res), function(Res) { data.frame(matrix(nrow = length(Africa_PAs_centroids),
  #     ncol = length(Var_URLs))) }) # List to hold a dataframe of EVI values per spatial resolution.
  #   for (Res in seq(length(Sites_res))) { # For each resolution of interest...
  #     for (Time in seq(length(Var_URLs))) { # For each raster, each containing EVI values for a time of a year...
  #       Var_df[[Res]][, Time] <- extract(subset(Var_data[[Res]], Time), Africa_PAs_centroids) # Extract EVI per site.
  #       print(paste('Res', as.character(Res), 'Time', as.character(Time), 'Complete')) }}; rm(Res, Time) # Progress.
  #   save(Var_df, file = 'EVI_df.rda') # Save "Var_df".
  #   load('EVI_df.rda') # Read in "Var_df".
  #   Var_df <- lapply(Var_df, function(Df) { return(Df[,-ncol(Df)]) }) # Remove last col of "Var_df", as refers to 2020.
  # 
  #   # To perform a null comparison to "Var_df", create a number of randomized versions of "Var_df", in which the EVI
  #   # values for each site are seasonally randomized, i.e., EVI values are randomly drawn from the distribution of all
  #   # values that refer to the same time of year (comment this out when the original "Var_df" will be analyzed):
  # 
  #   # Var_df_orig <- Var_df # Original copy of "Var_df" before it gets changed.
  #   # 
  #   # Seas_randomizer <- function(Vector, Indices, Method) { # For a vector of EVI data with seasonal indices...
  #   #   for (Index in unique(Indices)) { # For each unique seasonal index (e.g., index of EVI values of January 1-16)...
  #   #     Vals <- Vector[Indices == Index] # Identify the EVI values associated with the index.
  #   #     if (Method == 'Range') { Param <- range(Vals, na.rm = T); Func <- 'runif' } # Range of EVI at index. OR...
  #   #     if (Method == 'Norm') { Param <- c(mean(Vals, na.rm = T), sd(Vals, na.rm = T)); Func <- 'rnorm' } # Mean/sd.
  #   #     Vector[Indices == Index] <- do.call(Func, list(length(Vals), Param[1], Param[2])) } # Randomize EVI.
  #   #   return(Vector) } # Return the randomized vector of EVI values.
  #   # 
  #   # Null_number <- 100 # Number of random/null iterations to perform.
  #   # Null_out <- list() # List to hold the desired output from each null iteration.
  #   # for (Iter in 1:Null_number) { # For each iteration in which a randomized version of "Var_df" will be made...
  #   # 
  #   #   Var_df <- lapply(Var_df, function(Df) { as.data.frame(t(apply(Df, 1, Seas_randomizer, Indices = Var_periods,
  #   #     Method = 'Range'))) }) # Apply the above function to the rows of "Var_df" to make a null version of it.
  # 
  #   # Create a version of "Var_df" that refers to the scaled anomalies of each EVI value, to remove seasonal effects:
  #   
  #   Var_df_ano <- lapply(Var_df, function(Df) { # Per dataframe in "Var_df"...
  #     Df <- Df * 0.0001 # Scale the dataframe's EVI values according to the scale factor instructed by NASA.
  #     for (Row in 1:nrow(Df)) { # For each row/site in the dataframe...
  #       for (Index in unique(Var_periods)) { # For each unique seasonal index in "Var_periods"...
  #         Vals <- unlist(unname(Df[Row, Var_periods == Index])) # ID the site's EVI values associated with the index.
  #         Param <- c(mean(Vals, na.rm = TRUE), sd(Vals, na.rm = TRUE)) # Calculate the mean and sd of those EVIs.
  #         Df[Row, Var_periods == Index] <- (Vals - Param[1]) / Param[2] }} # Perform scaling using those calculations.
  #     return(Df) }) # Now that all dataframe values have been scaled to their seasonal equivalents, return the result.
  # 
  #   # Determine, from "Var_df_ano", a threshold to use for EVI values to be considered disturbances:
  # 
  #   Thr <- mean(sapply(Var_df_ano, function(Df) { # Per dataframe in "Var_df_ano"...
  #     quantile(as.vector(unlist(unname(Df))), 0.05, na.rm = TRUE) })) # Find its Xth percentile value. Then take mean.
  #   Thr <- -1 # Alternatively, manually set "Thr".
  # 
  #   # Plot a sample of rows in "Var_df_ano" as a scatterplot to determine how they are distributed across time:
  # 
  #   # par(mar = c(7,7,1,1)) # Set plot layout settings.
  #   # Rows <- sample(1:nrow(Var_df_ano[[1]]), 100) # Select rows in "Var_df_ano" to plot.
  #   # lapply(Var_df_ano, function(Df) { # Per dataframe in "Var_df_ano"...
  #   #   for (Row in 1:length(Rows)) { # Per row in the rows of the dataframe that will be plotted...
  #   #     if (Row == 1) { # If the first row is being plotted, create the baseline plot with axes...
  #   #       plot(unlist(unname(Df[Rows[Row],])), pch = 19, xlab = '', ylab = '', cex.axis = 2, cex = 0.3) } else {
  #   #       points(unlist(unname(Df[Rows[Row],])), pch = 19, cex = 0.3) }} # Else, add row on top of those prior.
  #   #   mtext('Time Index', side = 1, line = 4, cex = 3, font = 2) # Format the X axis label.
  #   #   mtext('EVI Anomaly', side = 2, line = 4, cex = 3, font = 2) # Format the Y axis label.
  #   #   abline(h = 0, col = 'firebrick2', lwd = 5) # Add horizontal line at Y = 0.
  #   #   abline(h = Thr, col = 'deepskyblue2', lwd = 5) }); rm(Rows) # Add horizontal line at Y = threshold.
  # 
  #   # Create a function that locates disturbance events for a given site in "Var_df_ano", where a disturbance is an EVI
  #   # anomaly value of less than a threshold that is isolated from other such values by a given span of time:
  # 
  #   Disturbance_finder <- function(Vals, Threshold, Span) { # For a vector of EVI values, a threshold, and a time span...
  #     if (sum(Vals < Threshold, na.rm = TRUE) > 0) { # If vector has EVI values below threshold (i.e., disturbances)...
  #       Disturbance_time <- which.min(Vals) # Find the time index of the highest-magnitude disturbance.
  #       Disturbance_mag <- Vals[Disturbance_time] # Find the magnitude of that disturbance.
  #       if (sum(Vals < Threshold, na.rm = TRUE) == 1) { # If only one disturbance (initially or after recursion)...
  #         return(c(Disturbance_time, Disturbance_mag)) } # Return that disturbance's time index and magnitude.
  #       else { # If more than one disturbance exists in "Vals", use recursion to process the next disturbance...
  #         Lower <- max(c(Disturbance_time-Span, 1)) # Determine the lower limit of EVI values surrounding disturbance.
  #         Upper <- min(c(Disturbance_time+Span, length(Vals))) # Determine the corresponding upper limit.
  #         Vals[Lower:Upper] <- NA # Remove the EVI values within the lower and upper limits of the present disturbance.
  #         return(c(Disturbance_time, Disturbance_mag, Disturbance_finder(Vals, Threshold, Span))) }}} # Recursion.
  # 
  #   # Create a function that determines the amount of time it takes for a given site to recover from a disturbance event:
  # 
  #   Recovery_timer <- function(Time, Row, Df, Counter = 1) { # For an index of time in a row of a "Var_df_ano" df...
  #     if (Time+4 > ncol(Df)) { return(NA) } # If the task below is not possible, return NA.
  #     Window <- (Time):(Time+3) # Determine the moving time window for which mean EVI values will be calculated.
  #     if (is.na(mean(unlist(unname(Df[Row,Window])), na.rm = TRUE))) { return(NA) } # Stop task here if not viable.
  #     if (mean(unlist(unname(Df[Row,Window])), na.rm = TRUE) > 0) { return(Counter) } # No. times until window mean +ve.
  #     else { return(Recovery_timer(Time+1, Row, Df, Counter+1)) }} # Continue recursively until window mean becomes +ve.
  # 
  #   # Use "Recovery_timer" to determine the "Span" argument to use in "Disturbance_finder". Specifically, determine the
  #   # span of time for which it takes X% of EVI anomalies that are independently a value of less than "Thr" to recover:
  # 
  #   Times <- c() # Vector to hold recovery time values of all anomalies.
  #   for (Df in Var_df_ano) { # For each dataframe in "Var_df_ano"...
  #     for (Site in 1:nrow(Df)) { # For each site in each dataframe...
  #       Disturbances <- which(Df[Site,] < Thr) # Determine which time points qualify as a disturbance.
  #       if (length(Disturbances) == 0) { next } # Skip this iteration of the loop if no disturbances exist.
  #       if (length(Disturbances) >= 2) { # If there are more than one disturbance...
  #         Unique <- sapply(seq(length(Disturbances)-1), function(D) { # Per disturbance that is not the last one...
  #           any(Df[Site, seq(Disturbances[D], Disturbances[D+1])] > 0, na.rm = TRUE) }) # Unique from one after it?
  #         Unique[length(Unique)+1] <- TRUE # The last disturbance is unique in any case, because none come after it.
  #         Disturbances <- Disturbances[Unique] } # Subset "Disturbances" down to those that are unique.
  #       if (length(Disturbances > 0)) { # If there is still at least one disturbance...
  #         Times <- c(Times, sapply(Disturbances, function(D) { Recovery_timer(D, Site, Df) })) } # Recovery time(s).
  #       print(paste0('Site ', Site, ' Complete')) }}; rm(Df, Site, Disturbances, Unique) # Record of progress.
  #   Sp <- unname(quantile(Times, 0.975, na.rm = TRUE)) # Xth percentile of "Times" values is the assigned "Sp" value.
  #   rm(Times) # This variable is no longer necessary and is large and computationally expensive to hold on to.
  # 
  #   # Per spatial resolution, calculate different metrics representing the ecosystem stability of each site, and record
  #   # the results of each metric in "Sites":
  # 
  #   for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
  # 
  #     # Determine the timings and magnitudes of EVI disturbances for each site, using "Disturbance_finder":
  # 
  #     Disturbances <- apply(Var_df_ano[[Dat]], 1, function(Row) { unname(Disturbance_finder(Row, Thr, Sp)) })
  # 
  #     # Calculate the resistance of each site, referring to its capacity to mitigate the occurrence of EVI disturbances:
  # 
  #     Sites[[Dat]][Var_name[grep('Resis', Var_name)]] <- sapply(1:length(Disturbances), function(Site) { # Per site...
  #       if (is.null(Disturbances[[Site]])) { return(NA) } else { # Only proceed if site has at least one disturbance.
  #         Resistances <- sapply(seq(2, length(Disturbances[[Site]]), by = 2), function(Val) { # Per disturbance val...
  #           1/abs(Disturbances[[Site]][Val]) }) # Calculate the resistance value associated with the disturbance.
  #         return(mean(Resistances, na.rm = TRUE)) } }) # Return the mean of the resistances across all disturbances.
  #     print(paste('Variable 1 complete for dataframe', Dat)) # Record of progress.
  # 
  #     # Calculate the recovery time of each site, referring to the time it takes to rebound from an EVI disturbance.
  #     # Concurrently, calculate the recovery rate of each site, referring to the rate of the disturbance rebound:
  # 
  #     Sites[[Dat]][Var_name[grep('Time|Rate', Var_name)]] <- t(sapply(1:length(Disturbances), function(Site) {
  #       if (is.null(Disturbances[[Site]])) { return(rep(NA, 2)) } else { # Proceed if site has at least one disturbance.
  #         Recoveries <- sapply(seq(1, length(Disturbances[[Site]]), by = 2), function(Val) { # Per disturbance time...
  #           Time <- Recovery_timer(Disturbances[[Site]][Val], Site, Var_df_ano[[Dat]]) # Calculate recovery time.
  #           Rate <- Disturbances[[Site]][Val+1] / Time # Calculate recovery rate.
  #           return(c(Time, Rate)) }) # Return the two calculations concatenated.
  #         return(rowMeans(Recoveries, na.rm = TRUE)) } })) # Return the means of calculations across all disturbances.
  #     print(paste('Variables 2 and 3 complete for dataframe', Dat)) # Record of progress.
  # 
  #     # Calculate the temporal invariability of each site, referring to its ability to maintain consistent levels of EVI:
  # 
  #     Sites[[Dat]][Var_name[grep('Invari', Var_name)]] <- sapply(1:nrow(Var_df[[Dat]]), function(Site) { # Per site...
  #       if (all(is.na(Var_df[[Dat]][Site, ]))) { return(NA) } else { # Only proceed if site has EVI values.
  #         Vals <- c() # Vector to hold the annual EVI values used to calculate the site's temporal stability.
  #         for (Year in seq(length(unique(Var_years)))) { # For each year in our time series...
  #           Tmp <- as.numeric(Var_df[[Dat]][Site, which(Var_years == unique(Var_years)[Year])]) # Obtain site's EVIs.
  #           Tmp <- Tmp * 0.0001 # Scale the EVI values according to the scale factor instructed by NASA.
  #           Tmp2 <- tryCatch({ PhenoDeriv(Tmp) }, # Determine start and end of season ("sos"/"eos") indices in EVI vals.
  #             error = function(E) { NA }) # Unless that is not possible (too many NA values), in which case return NA.
  #           if (all(is.na(Tmp2))) { Vals[Year] <- NA } else { # If no "sos"/"eos", leave annual EVI value as NA.
  #             if (Tmp2[names(Tmp2) == 'sos'] <= Tmp2[names(Tmp2) == 'eos']) { # If "sos" comes before "eos"...
  #               Vals[Year] <- max(Tmp[seq(Tmp2[names(Tmp2) == 'sos'], Tmp2[names(Tmp2) == 'eos'])], na.rm = T) } else {
  #               Vals[Year] <- max(Tmp[setdiff(1:length(Tmp), seq(Tmp2[names(Tmp2) == 'sos']-1, Tmp2[names(Tmp2) == 'eos']
  #                 +1))], na.rm = T) }} } # Get the maximum EVI value occurring in the "sos"-"eos" time period.
  #         Stability <- mean(Vals, na.rm = T) / sd(Vals, na.rm = T) # Calculate the site's overall temporal stability.
  #         return(Stability) } }) # Return the site's temporal stability value.
  #     print(paste('Variable 4 complete for dataframe', Dat)) }; rm(Dat, Disturbances, Thr, Sp) # Record of progress.
  # 
  #   # To analyze the null comparisons to "Var_df", complete the "Iter" "for" loop above and generate a histogram of the
  #   # outputs collected in "Null_out" (comment this out when the original "Var_df" will be analyzed instead):
  # 
  #   # Null_out[[Iter]] <- sapply(1:length(Sites), function(Idx) { cor(Sites[[Idx]][, Var_name[1]], abs(Sites[[Idx]][,
  #   #   Var_name[3]]), use = 'complete.obs', method = 'spearman') }) # Output of interest.
  #   # Var_df <- Var_df_orig # Reset "Var_df" to the original to prepare it for the next loop iteration.
  #   # print(paste('Iteration', as.character(Iter), 'complete')) }; rm(Iter) # Record of progress.
  #   # 
  #   # save(Null_out, file = 'EVI_null_method=range.rda') # Save "Null_out".
  #   # load('EVI_null_method=range.rda') # Read in "Null_out".
  #   # 
  #   # par(mar = c(4,8,1,1)) # Set plot layout settings.
  #   # sapply(1:length(Null_out[[1]]), function(Idx) { # For each set of vectors in "Null_out", plot a histogram...
  #   #   hist(sapply(Null_out, '[[', Idx), xlab = '', ylab = '', main = '', cex.axis = 2.5, col = 'cyan3', breaks = 15)
  #   #   mtext('Frequency', side = 2, line = 5, cex = 3.5) }) # Format the Y axis label.
  # 
  # # Collect data for all climate variables, which are also taken from a stack of raster layers:
  # 
  #   # Read in the necessary library for this procedure, and define the names of the variables:
  # 
  #   library(dismo) # For using the "biovars" function.
  #   Var_name <- c('Clim_TempMean',  'Clim_TempVar', 'Clim_PrecipTot', 'Clim_PrecipVar') # Variables of climate.
  # 
  #   # Process a series of URLs, where each URL points to a raster of a climatic condition, at 1-km resolution globally,
  #   # for a specific month in a specific year throughout the time series of interest. The URLs were acquired directly
  #   # via download from the CHELSA website from which the climate data is derived:
  # 
  #   Var_URLs <- lapply(list.files()[grep('Climate_URLs', list.files())], read.table) # URLs for min/max temp and precip.
  #   Var_URLs <- lapply(Var_URLs, function(List) { # Per URL list, each referring to a different climatic condition...
  #     Tmp <- List[apply(List, 1, function(Row) { grepl(paste(2000:2020, collapse = '|'), Row) }), ] # Subset to 2000+.
  #     Tmp <- as.character(Tmp) # Convert URLs to character.
  #     Tmp <- paste('/vsicurl/', Tmp, sep = '') # Add GDAL file system driver capability to each URL.
  #     Tmp <- Tmp[as.vector(sapply(1:(length(Tmp)/12), function(Year) { # Per year covered by the URLs...
  #       seq(Year, length(Tmp), length(Tmp)/12) }))] # Re-order URLs so that all those referring to the year are grouped.
  #     return(Tmp) }) # Return the processed URLs.
  # 
  #   # Per climatic condition, create a raster stack of all rasters represented by its URLs processed above. Then,
  #   # perform an analysis on each year in our time series, as follows: 1) obtain the stack of rasters, per condition,
  #   # representing that year; 2) crop each stack to Africa; 3) expand each to include a version of the coarser
  #   # resolution(s) of interest; and 4) for each site, extract and store its climate condition values from the stacks:
  # 
  #   Var_data <- lapply(Var_URLs, raster::stack) # Create a raster stack per climatic condition across all years.
  #   save(Var_data, file = 'Climate_stack.rda') # Save "Var_data", which contains a single raster stack per condition.
  #   load('Climate_stack.rda') # Read in "Var_data".
  #   
  #   Var_df <- lapply(1:length(Sites_res), function(Res) { # Per spatial resolution of interest...
  #     Df <- data.frame(matrix(nrow = length(Africa_PAs_centroids) * # Make a df with rows encompassing all sites...
  #       min(sapply(Var_URLs, length))/12, # ...as well as for all years for which there is climate data.
  #       ncol = length(Var_URLs)*12+1)) # Make columns representing each climatic condition for each month plus the year.
  #     colnames(Df) <- c(rep(c('MaxTemp', 'MinTemp', 'Precip'), each = 12), 'Year') # Set the df's column names.
  #     Df['Year'] <- rep(seq(2000, 2000-1+min(sapply(Var_URLs, length))/12), each = length(Africa_PAs_centroids))
  #     return(Df) }) # After setting the "Year" column above, return the finalized df.
  #   
  #   library(doParallel) # For setting up the parallel processing of the tasks below, to speed things up.
  #   Cluster <- makeCluster(detectCores() - 1) # Make a cluster with all available computer cores except one.
  #   registerDoParallel(Cluster) # Register the cluster so that it can be used.
  #   
  #   Process_raster <- function(Year) { # Create a function that performs the analysis.
  #     
  #     # Subset "Var_data" to include only the rasters representing the year number ("Year") of interest:
  #     
  #     Tmp <- lapply(1:length(Var_data), function(Condition) { subset(Var_data[[Condition]], 1:12 +
  #       12*(Year-1)) }) # "Year" is >= 1, so this discerns which indices of rasters in the stacks refer to the year.
  #     
  #     # Crop each raster stack in "Tmp" to the extent of Africa. Run this process in parallel to save time:
  #     
  #     Extent <- extent(Africa_PAs_rasters[[1]]) # The extent of Africa.
  #     Tmp1 <- foreach(Stack = iter(1:length(Tmp)), .packages = 'raster') %dopar% { # For each stack...
  #       crop(Tmp[[Stack]], Extent) } # Perform the cropping.
  #     print(paste0('Rasters for Year ', Year, ' have been cropped')) # Record of progress.
  #       
  #     # Per stack in "Tmp1", create a corresponding version of it in the coarser resolution(s) of interest:
  #       
  #     Tmp2 <- lapply(Tmp1, function(Cond) { aggregate(Cond, Sites_res[Sites_res > 1][1]) })
  #     print(paste0('Rasters for Year ', Year, ' have been aggregated')) # Record of progress.
  #     
  #     # Extract climate condition values for each site across "Tmp1" and "Tmp2", and store them in "Var_df":
  #     
  #     Values <- list(Tmp1, Tmp2) # Group all rasters with all values into a single list.
  #     for (Df in 1:length(Var_df)) { # Per df in "Var_df"...
  #       for (Cond in 1:length(Var_URLs)) { # Per climatic condition...
  #         Var_df[[Df]][1:length(Africa_PAs_centroids) + (length(Africa_PAs_centroids) * (Year-1)), 1:12 +
  #           (12 * (Cond-1))] <- extract(Values[[Df]][[Cond]], Africa_PAs_centroids) }} # Store in correct rows/cols.
  #     return(Var_df) } # Return the extracted values.
  #   
  #   for (Y in seq(1, min(sapply(Var_URLs, length))/12)) { # Per year with climate data...
  #     Var_df <- Process_raster(Y) }; rm(Y) # Apply the above function to that year.
  #   stopCluster(Cluster) # Shut down the cluster.
  #   
  #   save(Var_df, file = 'Climate_df.rda') # Save "Var_df".
  #   load('Climate_df.rda') # Read in "Var_df".
  # 
  #   # Convert the values in "Var_df" to degrees Celsius for temperature and kg/m^2 (mm) for precipitation:
  #   
  #   Var_df <- lapply(Var_df, function(Df) { # Per dataframe in "Var_df"...
  #     Df[, grep('Min|Max', colnames(Df))] <- Df[, grep('Min|Max', colnames(Df))]/10 - 273.15 # Temp unit conversion.
  #     Df[, grep('Prec', colnames(Df))] <- Df[, grep('Prec', colnames(Df))]/1000 # Precipitation unit conversion.
  #     return(Df) }) # Return the modified dataframe.
  #   
  #   # Use the climatic condition values processed above to derive climate variables and store them in "Sites":
  # 
  #   Var_df_biovars <- lapply(Var_df, function(Df) { # Per dataframe in "Var_df"...
  #     Df <- sapply(1:nrow(Df), function(Row) { # Per row in the dataframe...
  #       Tmp <- biovars(prec = unlist(unname(Df[Row, grep('Prec', colnames(Df))])), tmin = unlist(unname(Df[Row, grep(
  #         'Min', colnames(Df))])), tmax = unlist(unname(Df[Row, grep('Max', colnames(Df))]))) # Calculate bioclim vars.
  #       return(Tmp[colnames(Tmp) %in% c('bio1', 'bio4', 'bio12', 'bio15')]) }) # Keep only bioclim vars of interest.
  #     return(t(Df)) }) # Transpose "Df", with each row a site/year combination and each column a bioclim variable.
  #   
  #   Sites <- lapply(1:length(Sites), function(Idx) { # For each dataframe in "Sites"...
  #     Data <- Sites[[Idx]] # Create a variable equal to the dataframe.
  #     Data[Var_name] <- NA # Add in blank columns representing the derived climate variables.
  #     Data[Var_name] <- t(sapply(1:length(Africa_PAs_centroids), function(Site) { # Per site...
  #       Site_vals <- Var_df_biovars[[Idx]][seq(Site, nrow(Var_df_biovars[[Idx]]), by = length(Africa_PAs_centroids)), ]
  #       return(colMeans(Site_vals)) })) # Calculate the mean of its climate variables across all years.
  #     return(Data) }) # Return the dataframe, now with climate variables added to it.
  #   
  # # Collect data for all variables taken using an overlay of polygons (uncomment and process each variable one by one):
  # 
  # Var_name <- 'Env_Biome' # Variable to address each site's biome (forest, grassland, or desert).
  # Var_data <- shapefile('./official/wwf_terr_ecos.shp') # Shapefile for the variable.
  # for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
  #   Tmp <- over(Africa_PAs_centroids, Var_data) # Get the biome polygons that overlap with the sites.
  #   Sites[[Dat]][Var_name] <- ifelse(Tmp$BIOME %in% c(1:6, 12, 14), 'Forest', ifelse(Tmp$BIOME %in% 7:11, 'Grassland',
  #     'Desert')) }; rm(Dat, Tmp) # Record the biome as one of three options for each site, based on biome IDs.
  # 
  # library(exactextractr); library(sf) # For using the "exact_extract" and "st_as_sf" functions, respectively.
  # Var_name <- 'Mammal_Rich' # Variable to address each site's mammal species richness (and species identities).
  # Var_data <- shapefile('./MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp') # Shapefile for the variable.
  # Var_data <- Var_data[which(gIntersects(aggregate(Africa), Var_data, byid = TRUE)), ] # Subset to mammals in Africa.
  # Mammals <- list() # List to hold the names of mammals occurring at each site.
  # for (Res in seq(length(Sites))) { # For each dataframe/resolution in "Sites"...
  #   Tmp <- lapply(1:nrow(Var_data), function(Poly) { # For each mammal polygon...
  #     Df <- exact_extract(Africa_PAs_rasters[[Res]], st_as_sf(Var_data[Poly,]))[[1]] # Extract raster vals in polygon.
  #     Df <- Df[!is.na(Df[,1]), ] # Only maintain non-NA raster values.
  #     if (nrow(Df) != 0) { return(Df) } else { return(NULL) } }) # Return result, if has non-NA values.
  #   Mammals[[Res]] <- list() # List to hold mammal names for each site, for the resolution being addressed.
  #   for (Site in seq(length(Africa_PAs_centroids))) { # For each site...
  #     Tmp2 <- Africa_PAs_centroids@coords[Site, 2] # Reference the unique latitude of the site.
  #     Mammals[[Res]][[Site]] <- unique(Var_data@data$binomial[which(sapply(Tmp, function(Poly) { # For each polygon...
  #       if (is.data.frame(Poly)) { return(Tmp2 %in% Poly[which(Poly[,2] >= 0.05), 1]) } else { return(FALSE) } }))])
  #     } # Record if polygon overlaps with >=5% of site. Thus, identify the mammals occurring at the site.
  #   Sites[[Res]][Var_name] <- sapply(Mammals[[Res]], length) }; rm(Res, Tmp, Site, Tmp2) # Record mammal richness.
  # save(Mammals, file = 'Mammals.rda') # Save "Mammals".
  # 
  # # Collect data for all variables recorded at the country level (uncomment and process each variable one by one):
  # 
  # Var_name <- c('Hum_EPI', 'Hum_EPIchange') # Variables to address each site's country's EPI and EPI 10-year change.
  # Var_data <- read.csv('epi2020results20200604.csv') # CSV for the variables (EPI = Environmental Performance Index).
  # for (Dat in seq(length(Sites))) { # For each dataframe in "Sites", get EPI values per country...
  #   Sites[[Dat]][Var_name] <- Var_data[match(Sites[[Dat]]$PA_Country, Var_data$iso), c('EPI.new', 'EPI.change')]
  #   Sites[[Dat]][Sites[[Dat]]$PA_Country == 'SSD', Var_name] <- Sites[[Dat]][Sites[[Dat]]$PA_Country == 'SDN',
  #     Var_name] }; rm(Dat) # Because CSV doesn't recognize South Sudan, make its EPI values the same as those of Sudan.
  # 
  # Var_name <- 'Hum_GDP' # Variable to address each site's country's mean annual GDP per capita.
  # Var_data <- read.csv('API_NY/API_NY.GDP.PCAP.CD_DS2_en_csv_v2_2445354.csv', skip = 3) # CSV for the variable.
  # Tmp <- apply(Var_data[, grep(paste(as.character(2000:2019), collapse = '|'), colnames(Var_data))], 1, function(Row) {
  #   mean(Row, na.rm = TRUE) }) # Calculate the mean GDP per capita of each country from 2000-2019 (no 2020 data yet).
  # for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
  #   Sites[[Dat]][Var_name] <- Tmp[match(Sites[[Dat]]$PA_Country, Var_data$Country.Code)] }; rm(Dat, Tmp) # Input "Tmp".
  # 
  # save(Sites, file = 'Sites_V2_OneSigma+Full.rda') # Save "Sites" (could save the files as "Sites_V2_...", etc.).
  load('Sites_V2_OneSigma+Full.rda') # Read in "Sites".
  
# Per spatial resolution, collect data on each PA site that pertains to its component mammals:

  # Read in the data needed for this, including a dataset that was processed in Chapter 2 that contains synonyms of
  # the scientific names of mammal species:
    
  load('Mammals.rda') # Read in "Mammals", which contains lists of mammals occurring at each site.
  Mammals_synonyms <- read.csv('Mammals_NamesSynonyms_ProcessedFromPhylacine.csv') # Scientific name synonyms.
  Mammals_traits <- read.csv('Mammals_Traits_FromPhylacine.csv') # Mammal species' traits.
  
  # Update the scientific names in "Mammals" to be consistent with those in the Phylacine database:
  
  Mammals_namesnew <- list() # List to hold updated names of mammals.
  for (Res in 1:length(Mammals)) { # For each spatial resolution covered in "Mammals"...
    Mammals_namesnew[[Res]] <- list() # Make "Mammals_namesnew" a nested list.
    for (Site in 1:length(Mammals[[Res]])) { # For each site covered in "Mammals" within each resolution...
      Mammals_namesnew[[Res]][[Site]] <- gsub(' ', '_', Mammals[[Res]][[Site]]) # Re-format mammal names at the site.
      Mammals_namesodd <- Mammals_namesnew[[Res]][[Site]][!(Mammals_namesnew[[Res]][[Site]] %in%
        Mammals_synonyms$Binomial1)] # Find which, if any, mammal names at the site are inconsistent with Phylacine.
      if (length(Mammals_namesodd > 0)) { # If any mammal names at the site are inconsistent...
        for (Name in Mammals_namesodd) { # For each such name...
          Tmp <- which(Mammals_synonyms == Name, arr.ind = TRUE) # Which row of "Mammals_synonyms" has the name.
          if (nrow(Tmp) > 0) { # If the name is, in fact, in "Mammals_synonyms"...
            Mammals_namesnew[[Res]][[Site]][match(Name, Mammals_namesnew[[Res]][[Site]])] <- as.character(
              Mammals_synonyms$Binomial1[Tmp[1,1]]) # Update the name to match that used by Phylacine.
    }}}}}; rm(Res, Site, Mammals_namesodd, Tmp, Name)
  
  # Identify the mammal names that, after the above scientific name update, remain inconsistent with Phylacine:

  Mammals_namesano <- unique(unlist(sapply(1:length(Mammals_namesnew), function(Res) { # For each resolution...
    unique(unlist(Mammals_namesnew[[Res]]))[!(unique(unlist(Mammals_namesnew[[Res]])) %in% as.vector(as.matrix(
      Mammals_synonyms)))] }))) # Determine which names do not appear at all in "Mammals_syononyms".
  
  # Attempt to find the synonyms of the names in "Mammals_namesano" from taxonomic databases, synonyms that are consistent
  # with Phylacine (this code block was run once and commented out, with its outputs saved, to save time):
  
  # library(taxize) # For using the "synonyms" function.
  # Mammals_namesano_syn <- list() # List to hold synonyms of names in "Mammals_namesano" from online database search.
  # for (Name in 1:length(Mammals_namesano)) { # For each name...
  #   Mammals_namesano_syn[[Name]] <- tryCatch({ synonyms(Mammals_namesano[Name], db = 'itis') },
  #     error = function(X) {}) }; rm(Name) # Get name's synonyms from the "itis" db (other dbs didn't add anything).
  # 
  # Mammals_namesano_syn <- lapply(Mammals_namesano_syn, function(Name) { # Per name recognized by the db...
  #   if(!is.na(Name)) { return(c(names(Name), Name[[1]]$syn_name)) } else { return(NA) } }) # Keep name and synonyms.
  # Mammals_namesano_syn[sapply(Mammals_namesano_syn, function(Name) { is.character(Name) & length(Name) == 1 })] <-
  #   NA # Remove names that, despite being recognized by the db, are associated with no synonym.
  # 
  # Mammals_namesano_syn[[2]] <- c('Fukomys_mechowii', 'Fukomys_mechowi') # Correct this typing error.
  # Mammals_namesano_syn[grep('Colobus_caudatus', Mammals_namesano_syn)][[1]] <- Mammals_namesano_syn[grep(
  #   'Colobus_caudatus', Mammals_namesano_syn)][[1]][-2] # This specific synonym not in "Mammals_synonyms", so remove.
  # Mammals_namesano_syn[grep('Piliocolobus', Mammals_namesano_syn)][[1]][2] <- 'Piliocolobus_badius' # Make consistent.
  # Mammals_namesano_syn[grep('Nannospalax', Mammals_namesano_syn)][[1]][2] <- 'Spalax_ehrenbergi' # Make consistent.
  # Mammals_namesano_syn[grep('Acomys_selousi', Mammals_namesano_syn)] <- NA # Synonym not in "Mammals_synonyms".
  # Mammals_namesano_syn[grep('Otomys_karoensis', Mammals_namesano_syn)] <- NA # Synonym same as original, so remove.
  # Mammals_namesano_syn[grep('zanzibar', Mammals_namesano_syn)][[1]][2] <- 'Galagoides_zanzibaricus' # Make consistent.
  # 
  # Mammals_namesano_syn <- lapply(Mammals_namesano_syn, function(Name) { # Per remaining name with synonyms...
  #   if (!is.na(Name[[1]])) { # If the name exists...
  #     Name <- Name[1:2] # Only keep original name and its first synonym, the one that appears in "Mammals_synonyms".
  #     Name <- sub(' ', '_', Name) # Make consistent with "Mammals_synonyms" by replacing first space with "_".
  #     Name <- sub(' .*', '', Name) # Remove all sub-species from synonyms, as those are not in "Mammals_synonyms".
  #     return(Name) } else { return(Name) } }) # Return either processed name or original NA value.
  # 
  # save(Mammals_namesano_syn, file = 'Mammals_namesano_syn.rda') # Save "Mammals_namesano_syn".
  load('Mammals_namesano_syn.rda') # Read in "Mammals_namesano_syn".
  
  # Update various datasets now that the names of the mammals in our data have been updated as well as possible. Update
  # "Mammals_namesano" to contain only the species' names that at this point still lack a synonym that is consistent with
  # Phylacine. Update "Mammals_namesnew" to contain the appropriate synonyms of names. Finally, update the "Mammal_Rich"
  # columns in "Sites", as name updates may have changed the number of species occurring at certain sites:
  
  Mammals_namesano <- Mammals_namesano[sapply(Mammals_namesano_syn, function(Name) { length(Name) <= 1 })]
  Mammals_namesano_syn <- Mammals_namesano_syn[sapply(Mammals_namesano_syn, function(Name) { length(Name) > 1 })]
  
  for (Name in 1:length(Mammals_namesano_syn)) { # For each name to be updated in "Mammals_namesnew"...
    Mammals_namesnew <- rapply(Mammals_namesnew, function(List) { # For each list element in "Mammals_namesnew"... 
      if (Mammals_namesano_syn[[Name]][1] %in% List) { # If the list element contains the name to be updated...
        List[List == Mammals_namesano_syn[[Name]][1]] <- Mammals_namesano_syn[[Name]][2] # Update name.
        return(List) } else { return(List) }}, classes = 'character', how = 'list') }; rm(Name)
  
  Mammals_namesnew <- rapply(Mammals_namesnew, function(List) { unique(List) }, classes = 'character', how = 'list')
  for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
    Sites[[Dat]]['Mammal_Rich'] <- sapply(Mammals_namesnew[[Dat]], length) }; rm(Dat) # Update "Mammal_Rich" column.
  
  # Per spatial resolution, calculate community-level trait values for each site from each's component mammals (this code
  # block was run once and commented out, with its outputs saved in "Sites", to save time):
  
  # library(FD) # For using the "dbFD" function.
  # Mammals_FDnames <- c('FRic', 'FEve', 'FDiv', 'FDis', 'RaoQ') # Names of functional diversity metrics measured below.
  # 
  # for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
  #   for (Site in seq(nrow(Sites[[Dat]]))) { # For each site...
  #     if (Sites[[Dat]]$Mammal_Rich[Site] > 0) { # If the site has mammals...
  # 
  #       # Get the row indices of "Mammal_traits" that describe the mammals occurring at the site:
  # 
  #       Tmp <- match(Mammals_namesnew[[Dat]][[Site]], Mammals_traits$Binomial.1.2) # Rows of species at site.
  #       Tmp <- Tmp[!is.na(Tmp)] # Remove NA row indices.
  # 
  #       # Calculate basic community-level trait values, including the mean and SD of the body masses of the mammals at
  #       # the site, plus the means of the percentages of their diets that are composed of plants, verts, and inverts:
  # 
  #       Sites[[Dat]][Site, 'Mammal_MassMean'] <- mean(log(Mammals_traits$Mass.g[Tmp]/1000), na.rm = T) # Mean body mass.
  #       Sites[[Dat]][Site, 'Mammal_MassSd'] <- sd(log(Mammals_traits$Mass.g[Tmp]/1000), na.rm = T) # SD body mass.
  #       Sites[[Dat]][Site, 'Mammal_DietPlant'] <- mean(Mammals_traits$Diet.Plant[Tmp], na.rm = T) # Mean diet plant %.
  #       Sites[[Dat]][Site, 'Mammal_DietVert'] <- mean(Mammals_traits$Diet.Vertebrate[Tmp], na.rm = T) # Vert %.
  #       Sites[[Dat]][Site, 'Mammal_DietInvert'] <- mean(Mammals_traits$Diet.Invertebrate[Tmp], na.rm = T) # Invert %.
  # 
  #       # Calculate advanced community-level functional diversity metrics from the trait values of mammals at the site:
  # 
  #       Df <- data.frame(Trait1 = log(Mammals_traits$Mass.g[Tmp]/1000), Trait2 = Mammals_traits$Diet.Invertebrate[Tmp],
  #         Trait3 = Mammals_traits$Diet.Vertebrate[Tmp], Trait4 = Mammals_traits$Diet.Plant[Tmp]) # Mammal trait values.
  #       FuncDiv <- dbFD(Df, w.abun = F, stand.x = T, calc.CWM = F, messages = F) # Get metrics after scaling values.
  #       FuncDiv_names <- Mammals_FDnames[Mammals_FDnames %in% names(FuncDiv)] # Names of metrics.
  #       Sites[[Dat]][Site, paste('Mammal_', FuncDiv_names, sep = '')] <- unname(unlist(unname(FuncDiv[
  #         FuncDiv_names]))) } else { # Put metrics in "Sites", one column per metric.
  # 
  #       next }}}; rm(Dat, Site, Tmp, Df, FuncDiv, FuncDiv_names) # If the site does not have mammals, skip over it.

  # Per spatial resolution, calculate the phylogenetic diversity of mammals at each site (this code block was run once
  # and commented out, with its outputs saved in "Sites", to save time and space):
  
  set.seed(714) # For ensuring repeatability of all random results below.
  # library(phytools) # For using the "read.nexus" function, and general phylogenetic analyses.
  # 
  # Mammals_phylo <- read.nexus('Complete_phylogeny.nex') # 1,000 trees of mammals from Phylacine.
  # Mammals_phylo <- Mammals_phylo[sample(1:length(Mammals_phylo), 100, replace = FALSE)] # Randomly select 100 trees.
  # 
  # for (Dat in seq(length(Sites))) { # For each dataframe in "Sites"...
  #   for (Site in seq(nrow(Sites[[Dat]]))) { # For each site...
  #     Tmp <- vapply(Mammals_phylo, function(Tree) { # For each tree in "Mammals_phylo"...
  #       Pruned_tree <- keep.tip(Tree, Mammals_namesnew[[Dat]][[Site]][Mammals_namesnew[[Dat]][[Site]] %in%
  #         Tree$tip.label]) # Prune the tree down to only the mammals occurring at the site.
  #       return(sum(Pruned_tree$edge.length, na.rm = TRUE)) }, FUN.VALUE = numeric(1)) # Total branch lengths of tree.
  #     Sites[[Dat]][Site, 'Mammal_Phylo'] <- mean(Tmp, na.rm = TRUE) # Mean of branch length sums at site across trees.
  #     print(paste0('Dataframe ', Dat, ' Site ', Site, ' Complete')) }}; rm(Dat, Site, Tmp) # Record of loop's progress.

# Order the columns in the dataframes in "Sites" appropriately. Then save the dataframes as CSV files (based on downstream
# analysis choices, potentially save a subset of the dataframes):
  
for (Df in 1:length(Sites)) { # For each dataframe...
  Cols <- colnames(Sites[[Df]])[grep('Site', colnames(Sites[[Df]]))] # Columns not to be alphabetically ordered.
  Sites[[Df]] <- Sites[[Df]][, c(Cols, sort(setdiff(colnames(Sites[[Df]]), Cols)))] # Put "Cols" first, then order others.
  # write.csv(Sites[[Df]], paste0('Sites_', Sites_res[Df], 'km.csv'), row.names = FALSE)
  }; rm(Df, Cols)

###########################################################################################################################
###########################################################################################################################

# PART III: DATA MANIPULATION AND EXPLORATION
  
# For each relevant variable in "Sites", pair it with a full name to use for labelling purposes below:
  
Varnames <- list(c('Clim_PrecipTot', 'Annual Precipitation'), c('Clim_PrecipVar', 'Precipitation Seasonality'),
  c('Clim_TempMean', 'Mean Annual Temperature'), c('Clim_TempVar', 'Temperature Seasonality'), c('Eco_RecoveryRate',
  'Stability'), c('Eco_RecoveryTime', 'Resilience Inverse'), c('Eco_Resistance', 'Resistance'), c('Env_Biome', 'Biome'),
  c('Env_HetHab', 'Habitat Heterogeneity'), c('Env_HetTopo', 'Topographic Heterogeneity'), c('Hum_EPI',
  'Environmental Protection Index'), c('Hum_EPIchange', 'Change in EPI'), c('Hum_GDP', 'GDP Per Capita'),
  c('Hum_HII', 'Human Influence Index'), c('Hum_TravelTime', 'Travel Time to City'), c('Mammal_DietInvert',
  'Percentage of Diet Composed of Invertebrates'), c('Mammal_DietPlant', 'Percentage of Diet Composed of Plants'),
  c('Mammal_DietVert', 'Percentage of Diet Composed of Vertebrates'), c('Mammal_FDis', 'Mammal Functional Dispersion'),
  c('Mammal_FDiv', 'Mammal Functional Divergence'), c('Mammal_FEve', 'Mammal Functional Evenness'),
  c('Mammal_MassMean', 'Mammal Mean Body Mass'), c('Mammal_MassSd', 'Mammal Standard Deviation Body Mass'),
  c('Mammal_Phylo', 'Mammal Phylogenetic Diversity'), c('Mammal_RaoQ', 'Mammal Rao Q'), c('Mammal_Rich',
  'Mammal Species Richness'), c('PA_Area', 'Protected Area Size'), c('PA_Desig', 'Protected Area Designation'),
  c('PA_IUCN', 'Protected Area IUCN Category'), c('PA_Shape', 'Protected Area Shape'), c('Plant_LDMC',
  'Leaf Dry Matter Content'), c('Plant_LNC', 'Leaf Nitrogen Content'), c('Plant_LPC', 'Leaf Phosphorus Content'),
  c('Plant_SLA', 'Specific Leaf Area')) # Pairs of abbreviations and full names.
Varnames_abb <- sapply(Varnames, '[[', 1) # Abbreviated variable names.
Varnames_full <- sapply(Varnames, '[[', 2) # Full variable names.
Varnames_finder <- function(Names) { Varnames_full[match(Names, Varnames_abb)] } # Function to find full variable names.
  
# Choose the response variable of interest to use in downstream analyses. Each variable can be used and then compared:

Analysis_responsegroup <- 'Eco_' # Group of variables to be used as candidate response variables.
Analysis_responses <- colnames(Sites[[1]])[grep(Analysis_responsegroup, colnames(Sites[[1]]))] # Response variable options.
Analysis_response <- Analysis_responses[2] # Could be any number between 1 and length of "Analysis_responses".

# Visualize the correlations between all candidate response variables, before subsetting down to only the chosen one:

library(psych) # For using the "pairs.panels" function.
library(ggplot2) # For general plotting procedures.
par(mar = c(7,7,2,2)) # Adjust plot margins.

lapply(Sites, function(Df) { # Per dataframe in "Sites"...
  
  # Prepare a dataframe of all response variables of interest:
  
  Responses <- Df[, grep(Analysis_responsegroup, colnames(Df))] # Select columns referring to the response variables.
  Responses <- Responses[, -grep('Time|Invariability', colnames(Responses))] # Potentially subset "Responses".
  colnames(Responses) <- Varnames_finder(colnames(Responses)) # Modify column names.
  Responses <- Responses[complete.cases(Responses), ] # Keep only the rows of "Responses" that have no missing values.
  if (Analysis_responsegroup == 'Eco_') { # If the candidate response variables refer to ecosystem stability metrics...
    Responses <- sapply(Responses, function(Res) { abs(Res) }) } # Absolute-value each response variable.
  # Responses <- apply(Responses, 2, scale) # Z-score scale each column in "Responses".
  Responses <- Responses[, rev(sort(colnames(Responses)))] # Re-order the response variable columns.
  
  # Plot response variable correlations, doing so differently depending on the number of variables:
  
  if (ncol(Responses) > 2) { # If there are more than two response variables under consideration...
    
    pairs.panels(Responses, method = 'spearman', density = TRUE, ellipses = FALSE, hist.col = 'steelblue3', cex.cor = 5,
      cex.labels = 3, cex.axis = 2, pch = 19, gap = 0.4) } else { # Make a pair-wise correlation plot. OR, if only two...
  
    Responses <- as.data.frame(Responses) # Convert "Responses" to a data frame.
    print(cor(Responses, method = 'spearman')) # Print the correlation value of the response variables' relationship.
    
    Label_decs <- function(Lab) { sprintf('%.1f', Lab) } # Function to set axis tick labels to have one decimal digit.
    Label_breaks <- function(Resp) { seq(min(Resp), max(Resp), length.out = 5) } # Function to set number of tick labels.
    
    ggplot(Responses, aes_string(x = colnames(Responses)[1], y = colnames(Responses)[2])) + # Select data for plot.
      geom_point(size = 3) + geom_smooth(method = 'lm', col = 'deepskyblue3', lwd = 3) + # Scatterplot with best fit line.
      labs(x = colnames(Responses)[1], y = colnames(Responses)[2]) + # Set axis labels.
      scale_x_continuous(breaks = Label_breaks(Responses[,1]), labels = Label_decs) + # Format X axis tick labels.
      scale_y_continuous(breaks = Label_breaks(Responses[,2]), labels = Label_decs) + # Format Y axis tick labels.
      theme_classic() + # Make the plot have a classic-looking theme.
      theme(axis.text = element_text(size = 30, color = 'black'), axis.text.x = element_text(hjust = 0.5), axis.title =
        element_text(size = 35, color = 'black'), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0,
        l = 0)), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) } }) # Format axes.

# Visualize the spatial distribution of the chosen response variable's values in the dataframes in "Sites":

library(viridis) # For coloring the points plotted below.
library(fields) # For using the "image.plot" function.
library(terra) # For using the "sbar" function.
par(mfrow = c(1,1), mar = c(0,0,0,0)) # Set plot layout settings.

SiteResponse_mapper <- function(Df, Col_no, Map, Pts) { # Create a function that performs this visualization per dataframe.

  Vals <- abs(unlist(unname(Df[Analysis_response]))) # Obtain the response variable values from "Df".
  Cols <- rev(viridis(Col_no)) # Set the colors of the sites in the plot to be a color ramp of "Col_no" distinct colors.
  Breaks <- unname(quantile(Vals, seq(0, 1, length.out = Col_no + 1), na.rm = T)) # Breaks for map colors.
  for (Val in which(duplicated(Breaks))) { Breaks[Val] <- Breaks[Val-1] + 0.000001 }; rm(Val) # Prevent duplicate breaks.
  Cols_vals <- Cols[cut(Vals, Breaks, labels = F)] # Set the colors of sites based on their response variable values.

  plot(Map, col = 'gray80') # Plot Africa in the background.
  plot(Pts, col = Cols_vals, pch = 19, add = TRUE) # Plot the sites.
  image.plot(legend.only = TRUE, zlim = range(Vals), col = Cols, breaks = Breaks, axis.args = list(cex.axis = 2.5,
    yaxt = 'n'), legend.width = 0.5, legend.mar = 1.5, legend.shrink = 0.5) # Add legend with continuous colors to plot.
  sbar(2000, xy = c(-10,-40), lonlat = TRUE, adj = c(0.5, -0.6), lwd = 5, cex = 3) } # Scale bar.

lapply(Sites, function(Df) { SiteResponse_mapper(Df, 1000, Africa, Africa_PAs_centroids) }) # Apply the above function.

# Subset the dataframes in "Sites" down to the chosen response variable of interest:

library(dplyr) # For using the "relocate" function.
Sites_orig <- Sites # Save an original copy of "Sites" here, before manipulations to it are made below.
Sites <- lapply(Sites, function(Df) { # Per dataframe in "Sites"...
  Df <- Df[, -which(colnames(Df) %in% base::setdiff(Analysis_responses, Analysis_response))] # Remove unchosen responses.
  if (Analysis_responsegroup == 'Eco_') { # If the candidate response variables refer to ecosystem stability metrics...
    Df[, grep(Analysis_responsegroup, colnames(Df))] <- abs(Df[, grep(Analysis_responsegroup, colnames(Df))]) } # Abs val.
  Df <- Df %>% relocate(all_of(Analysis_response), .after = last_col()) # Place the response variable as the last column.
  return(Df) }) # Return the processed dataframe.
detach(package:dplyr) # Detach this package, as it may interfere with other functions used later on.

# Choose the mode of all downstream analyses, i.e., whether certain variables of interest will be included up-front in the
# first machine learning model made, or whether variables will be iteratively added into the model afterwards. Based on
# the mode of choice, choose key phrases that will be used in places below to select and remove variables, as appropriate:

Analysis_mode <- 'UpFront' # "UpFront" or "Iterative".
Analysis_itergroup <- ifelse(Analysis_responsegroup == 'Eco_', 'Mammal_', 'Eco_') # Group of variables to iterate.

if (Analysis_itergroup == 'Mammal_') { # If the group of iterative variables refer to mammal diversity metrics...
  Analysis_phrases <- ifelse(Analysis_mode == 'UpFront', 'FDis|Phylo', 'F|Phylo|RaoQ|Rich') } # Phrases for var selection.
if (Analysis_itergroup == 'Eco_') { # If the group of iterative variables refer to ecosystem stability metrics...
  Analysis_phrases <- ifelse(Analysis_mode == 'UpFront', 'Resis|Time|Inv', 'Resis|Recov|Inv') } # Phrases.

Analysis_eliminations <- ifelse(Analysis_mode == 'UpFront', paste('Site_|', Analysis_response, sep = ''), paste(
  Analysis_itergroup, '|Site_|', Analysis_response, sep = '')) # Phrases for variable removal in analyses below.
Analysis_nonpredictors <- paste('Site_|', Analysis_response, sep = '') # Phrases that refer to non-predictors of models.

# Choose whether all biomes will be analyzed together, or whether certain biomes will be focused on or removed:

Analysis_biome <- 'All' # "All", "Desert", "Forest", or "Grassland", depending on which biome(s) is/are of interest.
Analysis_biome_remove <- 'None' # "None", "Desert", "Forest", or "Grassland", depending on which biome(s) to remove.
Sites <- lapply(Sites, function(Df) { Df[, 'Env_Biome'] <- as.factor(Df[, 'Env_Biome']); return(Df) }) # Convert to factor.

# Based on downstream analyses, remove certain of the columns/variables and rows from each dataframe in "Sites":

Sites_removecols <- c(setdiff(grep('Hum_', colnames(Sites[[1]])), grep('GDP|HII|Travel', colnames(Sites[[1]]))), # Human.
  setdiff(grep('Plant_', colnames(Sites[[1]])), grep('LDMC$|SLA$', colnames(Sites[[1]]))), # Plant columns to remove.
  setdiff(grep(Analysis_itergroup, colnames(Sites[[1]])), grep(Analysis_phrases, colnames(Sites[[1]]))), # Iters to remove.
  grep('PA_IUCN|PA_Desig', colnames(Sites[[1]]))) # Any other random columns to remove.
Sites_removerows <- lapply(1:length(Sites), function(Df) { # Per dataframe in "Sites"...
  which(Sites_orig[[Df]][, 'Mammal_Rich'] <= 2 | !complete.cases(Sites[[Df]][, -Sites_removecols])) }) # ID rows to remove.

if (Analysis_biome != 'All') { # If a particular biome is to be analyzed in isolation, subset to the biome of interest...
  Sites_removerows <- lapply(Sites_removerows, function(Vec) { # Per vector of rows to remove...
    unique(c(Vec, which(Sites[[1]][, 'Env_Biome'] != Analysis_biome))) }) # Add rows not referring to biome of interest.
  Sites_removecols <- unique(c(Sites_removecols, grep('Biome', colnames(Sites[[1]])))) } # Remove biome column.
if (Analysis_biome_remove != 'None') { # If a biome is to be removed, remove that biome's rows...
  Sites_removerows <- lapply(Sites_removerows, function(Vec) { # Per vector of rows to remove...
    unique(c(Vec, which(Sites[[1]][, 'Env_Biome'] == Analysis_biome_remove))) }) } # Add rows of biome to be removed.

Sites <- lapply(1:length(Sites), function(Df) { return(Sites[[Df]][-Sites_removerows[[Df]], -Sites_removecols]) }) # Rm.

# Identify all variables in each dataframe in "Sites" that have negative values, and count the number of negative values
# in those that should have none. Based on literature/metadata, convert certain of those negative values as appropriate,
# or alternatively remove rows with negative values:

Sites_numvars <- colnames(Sites[[1]])[sapply(Sites[[1]], is.numeric)] # Vector of the names of all continuous variables.
Sites_catvars <- colnames(Sites[[1]])[sapply(Sites[[1]], is.factor)] # Vector of the names of all categorical variables.

Negvals_processor <- function(Df, Method) { # Create a function that performs this process per dataframe.
  Vars <- Sites_numvars[sapply(Sites_numvars, function(Var) { sum(Df[, Var] < 0) > 0 })] # Continuous vars with negatives.
  Vars <- Vars[grep('Plant', Vars)] # Of the variables in "Vars", only those plant-related should not be negative.
  print(sapply(Vars, function(Var) { sum(Df[, Var] < 0) })) # Print the number of negative values per variable in "Vars".
  if (Method == 'Convert') { # If negative values are to be converted, and not discarded...
    Df[, Vars] <- sapply(1:length(Vars), function(Idx) { # For each variable in "Vars"...
      Metadata <- Df[, Vars[Idx]] == floor(Df[, Vars[Idx]]) # Find the rows of the variable that are site metadata.
      Df[Metadata, Vars[Idx]] <- 0 # Convert those metadata values to 0.
      Df[, Vars[Idx]] <- abs(Df[, Vars[Idx]]) }) } # Convert the variable's remaining negative values to absolute value.
  else if (Method == 'Remove') { # If rows with negative values are to be removed entirely...
    Rows <- unname(which(apply(Df[, Vars], 1, function(Row) { any(Row <= 0) }))) # Identify rows with values <= 0.
    Df <- Df[-Rows, ] } # Remove those rows from the dataframe.
  return(Df) } # Return the dataframe following the negative value conversion or removal.

Sites <- lapply(Sites, Negvals_processor, Method = 'Remove') # Apply the above function to each dataframe in "Sites".

# Plot boxplots of all continuous variables in each dataframe in "Sites" to check for any major outliers. If outliers
# exist that should not (based on literature/metadata), remove the rows in each dataframe containing them:

Outlier_detector <- function(Df) { # Create a function that plots the boxplots per dataframe.
  Df <- scale(Df[, Sites_numvars]) # Z-score scale all continuous variables.
  Df <- Df[, -grep('Site', colnames(Df))] # Remove columns referring to site metadata.
  par(mar = c(15,5,0.5,0.5), mfrow = c(1,1)) # Set plot layout settings.
  boxplot(Df, xaxt = 'n', yaxt = 'n', pch = 20, col = 'darkolivegreen3', lwd = 1.5) # Make the boxplots.
  axis(1, at = 1:(ncol(Df)), labels = FALSE) # Format the X axis ticks.
  text(1:(ncol(Df)), par('usr')[3] - 0.6, labels = Varnames_finder(colnames(Df)), pos = 2, offset = 0, xpd = TRUE,
    srt = 90, cex = 1.25) # Format the X axis tick labels.
  axis(2, cex.axis = 2, font = 1) # Format the Y axis ticks.
  mtext('Scaled Value', side = 2, line = 3, cex = 2, font = 2) } # Format the Y axis label.

lapply(Sites, Outlier_detector) # Apply the above function to each dataframe in "Sites".

Outlier_rows <- lapply(Sites, function(Df) { # Per dataframe in "Sites"...
  Vars <- Sites_numvars[grep('HetHab|Plant_L', Sites_numvars)] # Identify variables with outliers to be removed.
  Rows <- unique(unlist(lapply(Vars, function(Var) { which(abs(Df[, Var]) * ifelse(grepl('HetHab', Var), 0.0001, 1) > 1)
    }))) }) # Per variable, identify its rows that contain outliers. Combine the results such that no row is repeated.

for (Vec in 1:length(Outlier_rows)) { # Per vector of rows in "Outlier_rows", each associated with a "Sites" dataframe...
  if (length(Outlier_rows[[Vec]]) > 0) { Sites[[Vec]] <- Sites[[Vec]][-Outlier_rows[[Vec]], ] }}; rm(Vec) # Remove rows.

# Visualize the distributions of all continuous variables in each dataframe in "Sites":

Distribution_visualizer <- function(Df, Type) { # Create a function that performs this visualization per dataframe.
  Cols <- Df[, Sites_numvars[-grep('Site_', Sites_numvars)]] # All numeric columns that are not site metadata.
  Plot_dim <- floor(sqrt(ncol(Cols))) # Determine the number of rows and columns that should be used to plot all plots.
  par(mfrow = c(Plot_dim, Plot_dim), mar = c(2,2,3,1)) # Set plot layout settings.
  Plots <- sapply(colnames(Cols), function(Col) { # Per column in "Cols"...
    if (Type == 'histogram') { # If the chosen type of plot is a histogram...
      hist(Cols[, Col], xlab = '', ylab = '', main = Varnames_finder(Col), cex.main = 1.2, col = 'cadetblue3') }
    else if (Type == 'density') { # If the chosen type of plot is a kernel density plot...
      plot(density(Cols[, Col], bw = 1), xlab = '', ylab = '', main = Varnames_finder(Col), cex.main = 1.2) # Density line.
      polygon(density(Cols[, Col], bw = 1), col = 'cadetblue3') } }) # Fill area under density line with color.
  return(Plots) } # Return the plots, now that they have been specified.

lapply(Sites, Distribution_visualizer, Type = 'histogram') # Apply the above function to each dataframe in "Sites".

# Visualize the relationships between all continuous variables in the dataframes in "Sites":

library(corrplot) # For using the "corrplot" function.

Relationship_visualizer <- function(Df) { # Create a function that performs this visualization per dataframe.
  Vars <- Sites_numvars[-grep('Site_', Sites_numvars)] # Remove "Site" variables from this analysis.
  colnames(Df) <- Varnames_finder(colnames(Df)); Vars <- Varnames_finder(Vars) # Modify variable names.
  Correlations <- cor(Df[, Vars], method = 'spearman', use = 'complete.obs') # Non-parametric correlations.
  par(mfrow = c(1,1)) # Set plot layout settings.
  return(corrplot(Correlations, type = 'lower', tl.col = 'black', tl.offset = 0.1, tl.cex = 0.75, cl.cex = 1)) }

lapply(Sites, Relationship_visualizer) # Apply the above function to each dataframe in "Sites".

# Visualize the relationship between a categorical variable of interest and the response variable in "Sites". Based on
# that, potentially transform the categorical variable by combining alike categories:

Sites_catvar <- Sites_catvars[-grep('Site', Sites_catvars)][1] # Choose categorical variable of interest.
if (Sites_catvar %in% colnames(Sites[[1]])) { # If the categorical variable of interest is present...
  par(mar = c(5,7,2,2), mfrow = c(1,1)) # Set plot layout settings.
  lapply(Sites, function(Df) { # Per dataframe in "Sites"...
    Response <- Df[, Analysis_response]; Predictor <- Df[, Sites_catvar] # Designate response and categorical variables.
    boxplot(Response ~ Predictor, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', pch = 20, lwd = 3, col = 'olivedrab3')
    axis(1, at = 1:(length(levels(Predictor))), labels = FALSE) # Format the X axis ticks.
    text(1:(length(levels(Predictor))), par('usr')[3], labels = levels(Predictor), pos = 1, offset = 1, xpd = TRUE,
      cex = 2.5) # Format the X axis tick labels.
    axis(2, cex.axis = 2, font = 1) # Format the Y axis ticks.
    mtext(Varnames_finder(Analysis_response), side = 2, line = 4, cex = 3, font = 2) }) } # Format the Y axis label.

# if ('Env_Biome' %in% colnames(Sites[[1]])) { # If the categorical variable of interest is present...
#   Sites <- lapply(Sites, function(Df) { # Per dataframe in "Sites", transform the categorical variable...
#     Df[, 'Env_Biome'] <- as.factor(ifelse(Df[, 'Env_Biome'] == 'Desert', 'Desert', 'NotDesert')) # Desert or not?
#     return(Df) }) } # Return the dataframe now that the transformation has been done.

# Perform variance inflation factor (VIF) analysis on each dataframe in "Sites" to assess collinearity between predictors:

library(car) # For using the "vif" function.

VIF_analyzer <- function(Df) { # Create a function that performs VIF analysis per dataframe.
  Predictors <- colnames(Df)[-grep(Analysis_eliminations, colnames(Df))] # List of predictors of interest.
  Formula <- paste(Analysis_response, '~', paste(Predictors, collapse = ' + ')) # Regression formula for VIF.
  return(vif(lm(Formula, data = Df))) } # Return the VIF value per predictor variable.

VIF <- lapply(Sites, VIF_analyzer) # Apply the above function to each dataframe in "Sites".

# Identify candidate interaction terms to be included in the machine learning models built below:

library(mgcv) # For using the "gam" function.

Interactions_numeric <- lapply(Sites, function(Df) { # Per dataframe in "Sites", find continuous interaction terms...
  Vars <- Sites_numvars[-grep(Analysis_nonpredictors, Sites_numvars)] # Remove metadata and response variable.
  Vars_paired <- combn(Vars, 2, simplify = FALSE) # Get every possible combination of variables in "Vars".
  Vars_significant <- sapply(1:length(Vars_paired), function(Pair) { # Per combination/pair of variables in "Vars"...
    Interaction <- Df[, Vars_paired[[Pair]][1]] * Df[, Vars_paired[[Pair]][2]] # Multiply the two to get their interaction.
    # Interaction <- as.vector(scale(Interaction)) # Scale the interaction term to standardize across pairs of variables.
    GAM <- gam(Df[, Analysis_response] ~ Interaction) # Run generalized additive model of the response vs. the interaction.
    Out_RSq <- unname(summary(GAM)$r.sq) # Extract the r-squared from the model.
    return(Out_RSq > 0.1) }) # Does the interaction explain a significant percentage of the variance in the response?
  return(Vars_paired[Vars_significant]) }) # Select interaction terms of significance.
Interactions_numeric <- Reduce(intersect, Interactions_numeric) # Choose terms selected for both dataframes in "Sites".

Interactions_numcat <- lapply(Sites, function(Df) { # Per dataframe in "Sites", find continuous-categorical terms...
  
  # Identify the continuous and categorical variables to be referenced below:
  
  Vars_num <- Sites_numvars[-grep(Analysis_nonpredictors, Sites_numvars)] # Relevant continuous variables.
  Vars_cat <- Sites_catvars[-grep(Analysis_nonpredictors, Sites_catvars)] # Relevant categorical variable(s).
  
  # For the categorical variable(s), determine if it exhibits a significant interaction with each continuous variable:
  
  if (length(Vars_cat >= 1)) { # If at least one categorical variable is to be analyzed...
    
    Var_cat <- Vars_cat[1] # Choose the categorical variable to be analyzed.
    Interactions <- list(); Counter <- 1 # List to store interactions, and loop counter.
    for (Var in Vars_num) { # For each relevant continuous variable...
      
      # Perform statistics to determine if the slope of the relationship of the continuous variable with the response
      # differs across the categorical variable's categories, a proxy for a categorical-continuous variable interaction:
      
      Stats <- summary(aov(as.formula(paste(Analysis_response, '~', Var, '*', Var_cat)), data = Df))
      Stats <- Stats[[1]]$`Pr(>F)`[3] # Extract p-value from statistical test.
  
      # In the case of the interaction being significant, visualize it using a scatterplot and record the interaction:
      
      if (Stats < 0.05) { # If the categorical-continuous interaction is significant...
        
        # Plot a scatterplot showing how the behavior of the continuous variable changes across categories (comment this
        # out if not done, as it may produce many plots):
        
        # print(ggplot(Df, aes_string(Var, Analysis_response, col = Var_cat)) + # Select data for plot.
        #   geom_point(size = 3) + geom_smooth(method = 'lm', size = 3) + # Set to be a scatterplot with a regression line.
        #   labs(x = Var, y = Analysis_response, col = gsub('(?!^)(?=[[:upper:]])', '\n', gsub('.*_', '', Var_cat),
        #     perl = TRUE)) + # Set axis and legend labels, using regex.
        #   scale_color_manual(labels = as.character(sort(unique(Df[, Var_cat]))), values = c('tan3', 'deepskyblue3',
        #     'plum4')[1:length(unique(Df[, Var_cat]))]) + # Set the scatter plot color scheme.
        #   theme_classic() + # Make the plot have a classic-looking theme.
        #   theme(axis.text = element_text(size = 30), axis.title = element_text(size = 35, face = 'bold'), legend.text =
        #     element_text(size = 30), legend.title = element_text(size = 35, face = 'bold'))) # Format axes and legend.
        
        # Record the interaction by storing the names of its variables as an element in "Interactions". Update the counter:
        
        Interactions[[Counter]] <- c(Var, Var_cat); Counter <- Counter + 1 }}
    
      return(Interactions) } else { return(NULL) } }) # Return interaction terms, if relevant.
Interactions_numcat <- Reduce(intersect, Interactions_numcat) # Choose terms selected for both dataframes in "Sites".

Interactions_total <- c(Interactions_numeric, Interactions_numcat) # Combine all candidate interaction terms.

# Determine the degree of spatial autocorrelation that exists in each continuous variable in each dataframe in "Sites"
# (comment this out if it is not to be done, as spatial autocorrelation corrections will be performed as needed anyways):

library(spatialRF) # For spatially-conscious machine learning analyses moving forward.

Distance_matrix <- pointDistance(Sites[[1]][, c('Site_Long', 'Site_Lat')], lonlat = TRUE) # Matrix of distances b/w sites.
Distance_matrix[upper.tri(Distance_matrix)] <- t(Distance_matrix)[upper.tri(Distance_matrix)] # Make matrix full.
Distance_matrix <- Distance_matrix/1000 # Convert distance matrix from meters to kilometers.

# SpatialAutocorrelation_examiner <- function(Df) { # Create a function that performs this analysis per dataframe.
#   Vars <- Sites_numvars[-grep('Site_', Sites_numvars)] # Remove latitude/longitude from the variables to be analyzed.
#   MoranI <- plot_training_df_moran(data = Df, predictor.variable.names = Vars[-length(Vars)], dependent.variable.name =
#     Analysis_response, distance.matrix = Distance_matrix) # Calculate Moran's I per continuous variable.
#   return(plot(MoranI)) } # Return a visual of the Moran's I outputs from above.
# 
# lapply(Sites, SpatialAutocorrelation_examiner) # Apply the above function to each dataframe in "Sites".

# Split each dataframe in "Sites" into training and testing subsets. Perform the split by stratifying it across space,
# such that the split is not spatially biased (i.e., no spatial data leakage exists between training and testing subsets).
# Before splitting, potentially Box-Cox transform the response variable in "Sites" to make its distribution more normal:

library(caret) # For general data manipulation and machine learning analyses moving forward.
library(blockCV) # For performing a spatially-conscious training-testing data split.

if (Analysis_response == 'Eco_Invariability') { # If the response variable refers to ecosystem invariability...
  Sites <- lapply(Sites, function(Df) { # Per dataframe in "Sites"...
    Transform_BoxCox <- preProcess(Df, method = list(BoxCox = Analysis_response)) # Set up the Box-Cox transformation.
    Df <- predict(Transform_BoxCox, newdata = Df) # Apply the Box-Cox transformation to the dataframe.
    return(Df) }) } # Return the dataframe now that the transformation has been completed.

TrainTest_splitter <- function(Df) { # Create a function that performs this split per dataframe.
  
  # Create a version of "Africa_PAs_centroids" that only contains the centroids/points pertaining to the rows in "Df":
  
  Centroids <- Africa_PAs_centroids[Africa_PAs@data$WDPA_PID %in% Df[, 'Site_ID']]
  
  # Perform the data split and store the results in a list:
  
  Split <- spatialBlock(Centroids, theRange = 300000, selection = 'systematic', showBlocks = F, progress = F, verbose = F)
  Split <- which(Split$foldID != 1) # Choose df rows for training.
  
  Train <- Df[Split, ]; Test <- Df[-Split, ] # Create the training and testing data subsets.
  List <- list(Split, Train, Test) # Create a list of the training subset row indices, as well as both subsets.
  names(List) <- c('Split', 'Train', 'Test') # Name the list elements for reference below.
  
  # Plot a map of Africa showing the sites designated for each of the training and testing subsets:
  
  par(mfrow = c(1,1), mar = c(0,0,0,0)) # Set plot layout settings.
  plot(Africa, col = 'gray80') # Plot Africa.
  Colors <- sapply(1:nrow(Df), function(Site) { ifelse(Site %in% Split, 'royalblue1', 'orange4') }) # Site colors.
  plot(Centroids, col = Colors, pch = 19, add = TRUE) # Plot the sites, colored by their designation.
  sbar(2000, xy = c(-20,-25), lonlat = TRUE, adj = c(0.5, -0.6), lwd = 5, cex = 3) # Plot scale bar.
  
  return(List) } # Return the finalized list made above.

Sites_split <- lapply(Sites, TrainTest_splitter) # Apply the above function to each dataframe in "Sites".

Sites_split_rows <- which(names(Sites_split[[1]]) == 'Split') # Index of "Sites_split" list elements with row indices.
Sites_split_train <- which(names(Sites_split[[1]]) == 'Train') # Index of "Sites_split" list elements with training subset.
Sites_split_test <- which(names(Sites_split[[1]]) == 'Test') # Index of "Sites_split" list elements with testing subset.

# Now that training and testing subsets have been taken from an un-transformed "Sites", transformations on the dataframes
# in "Sites" can now be performed. The same transformations will be performed on the training and testing subsets, but
# will be done within the machine learning process, to avoid data leakage between subsets and between data folds in the
# cross-validation process of model training. The potential transformations are as follows: Box-Cox transform continuous 
# predictors to make their distributions normal, and Z-score scale continuous predictors to put them on the same scale.
# The former is necessary to mitigate the effects of outliers, and the latter in case coefficient-based models are built:

Transformer <- function(Df) { # Create a function that performs this transformation per dataframe.
  
  # Devise the strategy by which columns/variables will be selected for Box-Cox transformation specifically:
  
  Strategy <- 'All' # "All" or "Selective". See below for details of each selection strategy.
  if (Strategy == 'All') { Cols <- Sites_numvars[-grep(Analysis_nonpredictors, Sites_numvars)] } # Choose all predictors.
  if (Strategy == 'Selective') { # Choose the predictors for transformation based on if they contain major outliers...
    Cols <- sapply(Sites_numvars, function(Col) { # Per continuous predictor column in the dataframe...
      !(grepl(Analysis_nonpredictors, Col)) & max(abs(scale(Df[, Col])), na.rm = TRUE) > 10 }) # Find if it has outliers.
    Cols <- names(Cols)[Cols] } # Get the names of the columns with outliers.
  
  # Perform the Box-Cox transformation on each variable selected:
  
  Transform_BoxCox <- preProcess(Df, method = list(BoxCox = Cols)) # Set up the transformation.
  Df <- predict(Transform_BoxCox, newdata = Df) # Apply it to the dataframe.
  
  # Perform the Z-score transformation on all non-metadata continuous predictor variables:
  
  # Cols_all <- Sites_numvars[-grep(Analysis_nonpredictors, Sites_numvars)]
  # Transform_ZScore <- preProcess(Df, method = list(center = Cols_all, scale = Cols_all)) # Set up the transformation.
  # Df <- predict(Transform_ZScore, newdata = Df) # Apply it to the dataframe.
  
  return(list(Df, Cols)) } # Return the transformed dataframe, as well as the variables that were Box-Cox transformed.

Transformer_output <- lapply(Sites, Transformer) # Apply the above function to each dataframe in "Sites".
Transformer_vars <- lapply(Transformer_output, '[[', 2)[[1]] # Extract the variables that were Box-Cox transformed.
Sites <- lapply(Transformer_output, '[[', 1) # Extract "Sites" from "Transformer_output".

# Generate spatial predictors for the dataframes in "Sites", as well as for each dataframe in the training and testing
# subsets in "Sites_split", in order to mitigate issues of spatial autocorrelation. These predictors are Moran's
# Eigenvector Maps (MEMs) derived from the distance matrix between sites (this code block was commented out, with its
# outputs saved, to save time):

SpatialPredictor_generator <- function(Split) { # Create a function that generates these predictors per set of sites.
  Distance_matrix_subset <- Distance_matrix[as.vector(Split), as.vector(Split)] # Subset matrix to sites of interest.
  MEM <- mem_multithreshold(distance.matrix = Distance_matrix_subset, distance.thresholds = round(seq(0, max(
    Distance_matrix_subset)/2, length.out = 3))) # Generate MEMs for different distance thresholds.
  MEM_rank <- rank_spatial_predictors(distance.matrix = Distance_matrix_subset, spatial.predictors.df = MEM) # Rank MEMs.
  MEM <- MEM[, MEM_rank$ranking] # Re-order the MEMs based on their ranking.
  return(MEM) } # Return the finalized MEMs.

# MEMs <- SpatialPredictor_generator(1:nrow(Distance_matrix)) # Apply the above function to all sites in "Sites".
# MEMs_traintest <- lapply(Sites_split, function(Element) { # Per list element in "Sites_split"...
#   Train_rows <- Element[[Sites_split_rows]] # Identify the rows of "Sites" used for training.
#   Test_rows <- setdiff(1:nrow(Distance_matrix), Train_rows) # Identify the rows of "Sites" used for testing.
#   lapply(list(Train_rows, Test_rows), function(Rows) { SpatialPredictor_generator(Rows) }) }) # Apply function.
# 
# save(MEMs, file = 'MEMs.rda'); save(MEMs_traintest, file = 'MEMs_traintest.rda') # Save "MEMs" and "MEMs_traintest".
load('MEMs.rda'); load('MEMs_traintest.rda') # Read in "MEMs" and "MEMs_traintest".

###########################################################################################################################
###########################################################################################################################

# PART IV: MACHINE LEARNING

# Establish the K-fold cross-validation mechanism that will possibly be inputted into each machine learning model. Ensure
# that the mechanism chooses the simplest model within one SD of the best model in order to avoid overfitting:

ML_CV <- trainControl(method = 'cv', number = 2, selectionFunction = 'oneSE', verboseIter = T, savePredictions = 'final')

# List out all predictors of interest that will be used in machine learning analysis:

ML_predictors <- sort(c(Sites_numvars[-grep(Analysis_eliminations, Sites_numvars)], # Numeric predictors of interest.
  colnames(Sites[[1]])[sapply(colnames(Sites[[1]]), function(Col) { # Per variable, determine which is/are...
    is.factor(Sites[[1]][, Col]) & !grepl('Site_', Col) })])) # ...categorical predictors that are of interest.
ML_predictors_orig <- ML_predictors # Save the original predictors in case some are removed in feature selection below.

# ML_predictors <- ML_predictors[!(grepl('TravelTime|HII|PA', ML_predictors))] # Perform general feature selection OR...
# ML_predictors <- c('Clim_PrecipTot', 'Env_HetHab', 'Mammal_FDis', 'Plant_LDMC') # ...more specific selection.
# Transformer_vars <- Transformer_vars[!(grepl('TravelTime|HII|PA', Transformer_vars))] # Same general selection OR...
# Transformer_vars <- ML_predictors # Same specific selection. 
ML_predictors_preprocess <- ML_predictors[!(grepl('Biome|IUCN|Desig', ML_predictors))] # Obtain pre-processable predictors.

ML_package <- 'caret' # "spatialRF" or "caret", depending on which package will be used to build models.
# if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...
#   ML_predictors <- sort(c(ML_predictors, sapply(Interactions_total, function(Int) { # Per chosen interaction(s)...
#     paste(Int[1], '*', Int[2]) }) )) } # Add it in the correct format.

# Run a linear model on every combination of predictors in "ML_predictors" and determine which combination leads to the
# lowest AICc value, representative of the best model fit. This will help to inform model formulas moving forward
# (comment this out if it is not to be done, as it takes time and may not be informative for non-linear models):

# library(MuMIn) # For using the "dredge" function below.
# ML_formula_options <- lapply(Sites, function(Df) { # Per dataframe in "Sites"...
#   dredge(lm(as.formula(paste(Analysis_response, '~', paste(ML_predictors, collapse = ' + '))), data = Df, na.action =
#     na.fail), rank = 'AICc') }) # Run linear model using every combo of "ML_predictors" to determine best predictor set.

# Create the baseline formula, referencing predictors and the response, that will be used in each machine learning model:

if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...
  ML_formula <- paste(Analysis_response, '~', paste(ML_predictors, collapse = ' + ')) } else { # Make equation. Or else...
  ML_formula <- c(ML_predictors, Analysis_response) } # List the names of the predictors and response.

# Train a model or model ensemble per training dataset, each representing a dataframe in "Sites". If necessary, add MEMs
# iteratively into the model formula until spatial autocorrelation is removed from the model's residuals. As part of
# model training, output partial dependence plots that show the effects of individual variables:

library(pdp) # For using the "partial" function to create partial dependence plots.
library(RRF); library(randomForest) # For the plain or regularized random forest model.

ML_model <- 'RRFglobal' # Name of model(s) recognized by the "caret" package.
ML_params <- expand.grid(mtry = 4, coefReg = 0.1) # Hyperparameters for a specific model.
ML_paramlen <- 2 # Alternative to model hyperparameters, the number of random hyperparameter combinations to try.

ML_varsapproach <- 'Select' # "All" or "Select" - should all or a subset of variables be visualized in ML outputs?
if (ML_varsapproach == 'All') { ML_varsinterest <- ML_predictors } else { # Choose either all variables or a subset...
  ML_varsinterest <- c('Clim_PrecipTot', 'Plant_LDMC', 'Env_HetHab', 'Mammal_FDis') }

ML_trainer <- function(Formula, Df, Model, Params, Len, Split, MEM_num, Res) { # Create a function that performs training.
  
  # Create a subsetted version of "Distance_matrix" that just includes the sites in the training data subset:
  
  Distance_matrix_train <- Distance_matrix[as.vector(Split), as.vector(Split)]
  
  # Using the training data subset, train the model using the arguments specified by the specific package being used:
  
  if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...
    BoxCox_vars <- Transformer_vars[!grepl(Analysis_eliminations, Transformer_vars)] # Variables to Box-Cox transform.
    Model_trained <- train(form = as.formula(Formula), data = Df, # Model formula and data source.
      preProcess = list(# center = ML_predictors_preprocess, scale = ML_predictors_preprocess, # Comment if not scaling.
        BoxCox = BoxCox_vars), # Box-Cox transform specified predictors.
      method = Model, trControl = ML_CV, tuneGrid = Params, # tuneLength = Len, # Model type and optimization settings.
      importance = TRUE) } else { # Train the model using caret, pre-processing the data in doing so. Or else...

    Df <- Transformer(Df)[[1]] # Apply the same data pre-processing techniques to "Df" as was done to "Sites" above.
    Model_trained <- rf_spatial(data = Df, dependent.variable.name = Formula[length(Formula)], predictor.variable.names =
      Formula[-length(Formula)], distance.matrix = Distance_matrix_train) # Train using data and distance matrix.
    Coords <- Df[, c('Site_Long', 'Site_Lat')]; colnames(Coords) <- c('x', 'y') # Site coords to be used in "rf_tuning".
    Model_trained <- rf_tuning(model = Model_trained, xy = Coords, repetitions = 10, num.trees = 750,
      mtry = 5, min.node.size = 10) } # Tune the model's hyperparameters to optimize its performance.
  
  # Plot partial dependence plots of each continuous variable included in the trained model to see how the values of the 
  # response change across the values of each variable, while holding the other variables constant:
  
  if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...

    Plots <- list(); Counter <- 1 # List to hold plots, and counter for the loop below.
    for (Predictor in ML_varsinterest[ML_varsinterest %in% Sites_numvars]) { # For each continuous variable...
      Plot_data <- partial(Model_trained, pred.var = Predictor) # Obtain the data for the partial dependence plot.
      if (Predictor == 'PA_Area') { Plot_data$PA_Area <- Plot_data$PA_Area / 1000 } # Specific correction for plotting.
      Plots[[Counter]] <- ggplot(Plot_data, aes_string(Predictor, 'yhat')) + # Begin the plot by specifying the data.
        geom_smooth() + # geom_point(size = 3) + # Set the plot to be a curve smoothed over points by LOESS regression.
        scale_y_continuous(breaks = seq(range(Plot_data$yhat)[1], range(Plot_data$yhat)[2], length.out = 10)[c(2,5,8)],
          labels = function(Label) { sprintf("%.3f", Label) }) + # Place and round y-axis tick values appropriately.
        scale_x_continuous(breaks = seq(range(Plot_data[, Predictor])[1], range(Plot_data[, Predictor])[2],
          length.out = 10)[c(2,5,8)], labels = function(Label) { sprintf("%.2f", Label) }) + # X-axis tick values.
        labs(x = paste(Varnames_finder(Predictor), ifelse(Predictor == 'PA_Area', '(/1000)', ''))) + # X label.
        theme_classic() + # Make the plot have a pleasant-looking theme.
        theme(axis.text = element_text(size = ifelse(ML_varsapproach == 'All', 10, 15), color = 'black'), axis.title.x =
          element_text(size = ifelse(ML_varsapproach == 'All', 12, 17.5), margin = ggplot2::margin(t = 10, r = 0,
          b = 0, l = 0)), axis.title.y = element_blank()) # Format axes.
      print(paste0('Plot ', Counter, ' of ', length(ML_varsinterest[ML_varsinterest %in% Sites_numvars]), ' Done'))
      Counter <- Counter + 1 } # Update the loop counter.

    if (ML_varsapproach == 'All') { do.call('grid.arrange', Plots) } else { # Arrange all var plots on same panel. OR...
      do.call('grid.arrange', c(grobs = lapply(Plots, '+', theme(plot.margin = ggplot2::margin(10,250,10,250))),
        ncol = 1)) }} # Arrange the plots of those variables that were selected to be of interest.
  
  # Calculate Moran's I of the trained model's residuals to assess if spatial autocorrelation exists across them:
  
  if (ML_package == 'caret') { Model_res <- residuals(Model_trained) } else { Model_res <- Model_trained$residuals$values }
  Model_moran <- moran(x = Model_res, distance.matrix = Distance_matrix_train, verbose = FALSE)
  
  # If spatial autocorrelation in the model's residuals didn't exist at all, has been dealt with manually using the
  # "caret" package, or has been handled automatically using "spatialRF", return the trained model, its residuals' Moran's
  # I, and the number of MEMs that were potentially added to the model. Otherwise, manually handle the autocorrelation
  # by iteratively adding MEMs into the model, and only then return the outputs of interest:
  
  if (Model_moran$test$p.value >= 0.00 | ML_package == 'spatialRF') { # If autocorrelation is no longer an issue...
    return(list(Model_trained, Model_moran$test, ncol(MEMs) - MEM_num)) } else { # Return outputs. Or else...
    MEMs_train <- MEMs_traintest[[Res]][[1]] # Reference the MEMs that are associated with the training dataset.
    Df[, colnames(rev(MEMs_train))[MEM_num]] <- rev(MEMs_train)[, MEM_num] # Add MEM (highest-ranked first) to df.
    Formula <- paste(Formula, '+', colnames(rev(MEMs_train))[MEM_num]) # Add that MEM to the model formula.
    return(ML_trainer(Formula, Df, Model, Params, Len, Split, MEM_num-1, Res)) }} # Recursively apply function on next MEM.

Sites_trained <- lapply(1:length(Sites_split), function(Res) { # For each resolution of "Sites" being analyzed...
  return(ML_trainer(ML_formula, # Apply the above function, the first argument being the model formula, then...
    Sites_split[[Res]][[Sites_split_train]], # ...the training data subset...
    ML_model, ML_params, ML_paramlen, # ...the model of choice and associated model parameters OR parameter combinations...
    Sites_split[[Res]][[Sites_split_rows]], # ...the row indices of the training data subset...
    ncol(MEMs), Res)) }) # ...the number of MEMs potentially included in the model, and the resolution being analyzed.

# Use the trained model per training dataset to make predictions on the corresponding testing datasets:

Sites_trained_model <- which(sapply(Sites_trained[[1]], function(Item) { # Per model trained above...
  TRUE %in% grepl('train|rf', class(Item)) })) # Find the index within "Sites_trained" of the actual trained model itself.

ML_tester <- function(Model, Df) { # Create a function that performs model predictions for a dataframe of interest.
  
  # Obtain the predictions. For the "spatialRF" package, ensure that data preprocessing is done as it was in training:
  
  if (ML_package == 'spatialRF') { Df <- Transformer(Df)[[1]] } # Data pre-processing - same as training subset above.
  Preds <- predict(Model, Df) # Use the chosen model to make predictions on the dataframe of interest.
  if (ML_package == 'spatialRF') { Preds <- Preds$predictions } # Reference predicted values from testing.
  
  # Post-process predictions to improve their robustness. Do this by determining the features of the sites for which
  # bad predictions were made, and then correcting for those bad predictions:
  
  # if (Analysis_response == 'Eco_RecoveryRate') { # For the recovery rate response variable...
  #   Preds[Df[, 'Clim_PrecipTot'] < 150] <- Preds[Df[, 'Clim_PrecipTot'] < 150] * 0.8 } # Low-precip sites overpredicted.
  # if (Analysis_response == 'Eco_Resistance') { # For the resistance response variable...
  #   Preds[Df[, 'Env_Biome'] == 'Desert' & Df[, 'Mammal_FDis'] > 1.9] <- Preds[Df[, 'Env_Biome'] == 'Desert' &
  #     Df[, 'Mammal_FDis'] > 1.9] * 1.1 } # High-FD sites that are deserts are underpredicted.

  # Return the final, unnamed vector of predictions of the response variable made on the testing dataset:
  
  return(unname(Preds)) }

Sites_tested <- lapply(1:length(Sites_trained), function(Res) { # For each resolution being analyzed...
  ML_tester(Sites_trained[[Res]][[Sites_trained_model]], Sites_split[[Res]][[Sites_split_test]]) }) # Predict.

# Evaluate each model per resolution represented by "Sites" by plotting the actual vs. predicted response variable values:

ML_evaluator <- function(Df, Rows, Predicted, Actual) { # Create a function that performs model evaluation per resolution.
  
  # Obtain evaluation metrics like mean squared error and r-squared:

  Evaluate <- postResample(pred = Predicted, obs = Actual)
  Evaluate['RMSE_Norm'] <- abs(Evaluate['RMSE'] / mean(Actual, na.rm = TRUE))

  # Plot the actual vs. predicted response variable values as a scatterplot with a regression line and perfect-model line:
  
  par(mfrow = c(1,1), mar = c(8,8,2,2)) # Set plot layout settings.
  plot(Predicted, Actual, xlab = '', ylab = '', xlim = range(c(Actual, Predicted), na.rm = TRUE), ylim = range(c(Actual,
    Predicted), na.rm = TRUE), pch = 19, cex.axis = 2.5, cex = 3) # Make the plot, equalizing the ranges of the axes.
  mtext(paste('Predicted', Varnames_finder(Analysis_response)), side = 1, line = 5, cex = 3) # X axis label.
  mtext(paste('Actual', Varnames_finder(Analysis_response)), side = 2, line = 5, cex = 3) # Y axis label.
  
  lines(x = c(-10000,10000), y = c(-10000,10000), col = 'deepskyblue3', lwd = 10, lty = 'dashed') # Perfect-model line.
  
  Regression <- lm(Actual ~ Predicted) # Perform linear regression on the above plot.
  abline(Regression, col = 'darkgoldenrod3', lwd = 10) # Plot the resultant regression line on the plot.
  
  mtext(paste('Y = ', as.character(round(Regression$coefficients[2], 2)), 'X + ', as.character(round(
    Regression$coefficients[1], 2)), sep = ''), side = 1, line = -31, cex = 2, font = 1, at = range(c(Actual,
    Predicted))[1], adj = 0.05) # Equation with slope and intercept in plot.
  mtext(paste('Norm RMSE =', as.character(round(Evaluate['RMSE_Norm'], 2))), side = 1, line = -28, cex = 2, font = 1,
    at = range(c(Actual, Predicted))[1], adj = 0.05) # Normalized RMSE value in plot.
  mtext(paste('R-Squared =', as.character(round(Evaluate['Rsquared'], 2))), side = 1, line = -25, cex = 2, font = 1,
    at = range(c(Actual, Predicted))[1], adj = 0.05) # R-squared value in plot.
  
  return(Evaluate) } # Return the evaluation metrics above.

Sites_evaluated <- lapply(1:length(Sites_tested), function(Res) { # For each resolution of "Sites" being analyzed...
  return(ML_evaluator( # Apply the above function, the arguments being...
    Sites[[Res]], # The dataframe addressing the resolution of interest...
    Sites_split[[Res]][[Sites_split_rows]], # ...the row indices of the training data subset...
    Sites_tested[[Res]], # ...the values of the response variable predicted by the model...
    Sites_split[[Res]][[Sites_split_test]][, Analysis_response])) }) # ...and the corresponding actual values.

# Visualize the importance of the variables that went into the optimized model per resolution represented by "Sites". In
# so doing, determine if each variable had a statistically significant impact on the model's performance, either by
# making boxplots of multiple simulations of variable importance, or by a self-made permutation test:

library(gtools) # For using the "permute" function.
ML_importancer <- function(Variable, Df, Res, Metric) { # Create a function that statistically evaluates each variable.
  Df[, Variable] <- permute(Df[, Variable]) # Permute the variable of interest by jumbling the order of its values.
  Tested <- ML_tester(Sites_trained[[Res]][[Sites_trained_model]], Df) # Make predictions on "Df" after permutation.
  Evaluated <- postResample(pred = Tested, obs = Df[, Analysis_response]) # Evaluate those predictions.
  return(unname(Evaluated[Metric])) } # Return the r-squared or other metric value of the evaluations.

library(vip) # For using the "vi_permute" function.
library(reshape) # For using the "melt" function.
ML_simnum <- 10 # Number of simulations to perform in the function below.

lapply(1:length(Sites_trained), function(Res) { # For each resolution of "Sites" being analyzed...
  
  # Select, process, and analyze the variable importance values and labels that will be plotted next:
  
  Idx <- Sites_trained_model # Trained model list index.
  
  if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...
    
    # Obtain a dataframe of variable importance values. Do so by permuting each variable "nsim" number of times, and
    # determining, per simulation, the degree to which a model evaluation metric changes because of the permutation.
    # Remove all interaction terms from this procedure (at least for now):
    
    Data <- vi_permute(Sites_trained[[Res]][[Idx]], train = Sites_split[[Res]][[Sites_split_train]], feature_names =
      ML_predictors[!grepl('\\*', ML_predictors)], target = Analysis_response, metric = 'rsquared', pred_wrapper =
      predict, nsim = ML_simnum) # Perform the permutation per variable.
    Data <- attr(Data, 'raw_scores') # Get the raw importance values from the simulations above.
    Data <- melt(Data) # Reshape "Data" to have a column representing variables and another with their importance values.
    Data <- Data[, c('X1', 'value')] } # Subset "Data" to its appropriate columns.
      
  if (ML_package == 'spatialRF') { # If the package being used to build machine learning models is "spatialRF"...
    Data <- Sites_trained[[Res]][[Idx]]$importance$per.variable } # Obtain a data frame of variable importance values.
  
  colnames(Data) <- c('Variable', 'Importance') # Make the column names of "Data" cleaner and compatible.
  
  Data$Group <- ifelse(grepl('Clim_', Data$Variable), 'Climate', ifelse(grepl('Env_', Data$Variable), 'Environment',
    ifelse(grepl('Eco_', Data$Variable), 'Ecosystem Stability', ifelse(grepl('Plant_', Data$Variable),
    'Plant Diversity', ifelse(grepl('Hum_', Data$Variable), 'Human Impacts', ifelse(grepl('PA_', Data$Variable),
    'Protected Area', 'Mammal Diversity')))))) # Assign labels to groups.
  
  Colors <- data.frame(Group = sort(unique(Data$Group))) # Initiate a data frame to hold plot colors by group.
  Colors$Col <- ifelse(grepl('Climate', Colors$Group), 'cornflowerblue', ifelse(grepl('Environment', Colors$Group),
    'gold3', ifelse(grepl('Ecosystem', Colors$Group), 'peachpuff2', ifelse(grepl('Plant', Colors$Group),
    'olivedrab', ifelse(grepl('Human', Colors$Group), 'orange3', ifelse(grepl('Protected', Colors$Group), 'slategray4',
    'mediumpurple3')))))) # Assign colors to groups.
  
  if (ML_package == 'spatialRF') { # If the package being used to build machine learning models is "spatialRF"...
    Data_met <- 'Rsquared' # Metric to use for the significance tests performed below.
    Data$Significant <- sapply(as.character(Data$Variable), function(Var) { # Per variable to be evaluated...
      if (grepl('Biome', Var)) { Var <- 'Env_Biome' } # Resolve categorical variable naming issue. Make null distribution
        # next. (NOTE: the above line of code may need to be expanded to include other categorical variables, as needed).
      Null_values <- replicate(ML_simnum, ML_importancer(Var, Sites_split[[Res]][[Sites_split_test]], Res, Data_met))
      return(Sites_evaluated[[Res]][Data_met] - quantile(Null_values, 0.95) > 0 & # Is test statistic in the tail AND...
        Sites_evaluated[[Res]][Data_met] - mean(Null_values, na.rm = T) > 0.0001) }) } # ...sufficiently better than null?
  
  if (ML_package == 'caret') { Data$Points <- 0 }

  Data$Name <- Varnames_finder(Data$Variable) # Get the full name of each variable in the "Variable" column.
  Data <- Data[Data$Variable %in% ML_varsinterest, ] # Subset "Data" to the variables of interest.
  # Data$Variable <- gsub(':(.*)_', ':', Data$Variable) # Remove group codes from interaction terms in "Variable" column.
  # Data$Variable <- gsub('.*_', '', Data$Variable) # Remove remaining group codes from all terms in "Variable" column.
  
  Ranks <- data.frame(Name = sort(unique(Data$Name))) # Initiate a data frame to hold importance ranks per group.
  Ranks$Median <- sapply(Ranks$Name, function(N) { median(Data$Importance[Data$Name == N], na.rm = TRUE) }) # Group meds.
  Ranks <- Ranks[rev(order(Ranks$Median)), ] # Re-order data frame by median values.
  Ranks$Group <- sapply(Ranks$Name, function(N) { Data$Group[Data$Name == N][1] }) # Pair variables with their groups.
  print(lapply(unique(Ranks$Group), function(G) { list(G, mean(which(Ranks$Group == G))) })) # Print mean rank per group.
  
  # Plot the variable importance values, as well as their labels and colors, as a bar graph or boxplot:
  
  if (ML_package == 'caret') { # If the package being used to build machine learning models is "caret"...
    print(ggplot(Data, aes(x = reorder(Name, Importance), y = Importance, color = Group)) + # Get data.
      geom_boxplot(fill = 'gray90', lwd = 1, fatten = 1) + # Set the plot to be a boxplot.
      scale_color_manual(values = Colors$Col) + # Manually select the colors to associate with each group in the plot.
      labs(y = paste('Importance for', Varnames_finder(Analysis_response))) + # Set plot's axis labels.
      theme_minimal() + # Make the plot have a minimal-looking theme.
      theme(axis.text.x = element_text(size = 25, color = 'black'), axis.text.y = element_text(size = 18, color = 'black'),
        axis.title.x = element_text(size = 30, margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)), axis.title.y =
        element_blank(), legend.text = element_text(size = 18, color = 'black'), legend.title = element_text(size = 20)) +
      coord_flip()) } # Flip the plot such that the boxes point horizontally, instead of vertically.
  
  if (ML_package == 'spatialRF') { # If the package being used to build machine learning models is "spatialRF"...
    print(ggplot(Data, aes(x = reorder(Name, Importance), y = Importance, fill = Group, shape = Significant)) + # Data.
      geom_bar(stat = 'identity', color = 'black') + # Set the plot to be a bar plot with select colors.
      scale_fill_manual(values = Colors$Col) + # Manually select the colors to associate with each group in the plot.
      geom_point(show.legend = FALSE, size = 2, stroke = 2) + # Put asterisks next to statistically significant variables.
      scale_shape_manual(values = c(NA, 8)) + # Ensure the points are asterisks and only appear for the right variables.
      labs(x = 'Attribute', y = paste('Importance for', Varnames_finder(Analysis_response))) + # Set plot's axis labels.
      theme_minimal() + # Make the plot have a minimal-looking theme.
      theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 18), axis.title = element_text(
        size = 30, face = 'bold'), legend.text = element_text(size = 18), legend.title = element_text(size = 20,
        face = 'bold')) + # Format axes and legend.
      coord_flip()) } }) # Flip the plot such that the bars point horizontally, instead of vertically.

# Use the above modeling framework to make an ensemble of models and compare their performances to each other:

ML_comparer <- function(Res, Term) { # Create a function that builds models, each with a term added to the formula above.
  
  # Perform the modeling steps as was done above, using the functions that were created above to do so:
  
  if (ML_package == 'caret') { Comp_formula <- paste(ML_formula, '+', Term) } else { # Add term to "ML_formula" eqn. Or...
    Comp_formula <- sort(c(ML_formula, Term)) } # ...add the term to the list of variables in "ML_formula".
  Comp_trained <- ML_trainer(Comp_formula, Sites_split[[Res]][[Sites_split_train]], ML_model, ML_params, ML_paramlen,
    Sites_split[[Res]][[Sites_split_rows]], ncol(MEMs), Res) # Train model.
  Comp_tested <- ML_tester(Comp_trained[[Sites_trained_model]], Sites_split[[Res]][[Sites_split_test]]) # Test model.
  Comp_evaluated <- postResample(pred = Comp_tested, obs = Sites_split[[Res]][[Sites_split_test]][, Analysis_response])
  
  # Create, process, and return a sequence of outputs of interest from the modelling steps performed above:
  
  Comp_out <- c(Comp_trained[[which(sapply(Comp_trained, is.data.frame))]]$p.value, unname(Comp_evaluated['Rsquared']))
  names(Comp_out) <- c('p-value', 'R-Squared') # Label the key outputs obtained in the line above.
  return(Comp_out) } # With the key outputs obtained and labeled, return them.

ML_predictors_new <- setdiff(Sites_numvars, ML_predictors_orig) # Identify predictors to add as terms to the model.
ML_predictors_new <- ML_predictors_new[!(grepl(Analysis_nonpredictors, ML_predictors_new))] # Non-metadata predictors.

if (length(ML_predictors_new) > 0) { # If such predictors of interest exist...
  Sites_modcompared <- lapply(1:length(Sites), function(Res) { # For each resolution of "Sites" being analyzed...
    lapply(ML_predictors_new, function(Term) { ML_comparer(Res, Term) }) }) } # Make model with each new term added in.

# Plot the results obtained from "Sites_modcompared" to visualize comparisons in the ensemble of models' performances:

ML_compareplotter <- function(Idx) { # Create a function that makes this plot per resolution of "Sites" being analyzed.
  
  # Compile the data to be visualized, namely the R-squared of each model made above and the labels of those models:
  
  Data <- sapply(Sites_modcompared[[Idx]], '[[', 2) # Extract R-squared values from each element of "Sites_modcompared".
  names(Data) <- ML_predictors_new # Associate each R-squared value with the term that was added into its model.
  Data <- c(Sites_evaluated[[Idx]][[grep('squared', names(Sites_evaluated[[Idx]]))]], Data) # Add in baseline model R-sq.
  names(Data)[1] <- 'Control' # Label the baseline model's R-squared value.
  
  # Make a bar plot of each R-squared value listed in "Data":
  
  par(mfrow = c(1,1), mar = c(11,11,2,2)) # Set plot layout settings.
  Bar <- barplot(Data, xaxt = 'n', cex.axis = 3, col = ifelse(grepl('Control', names(Data)), 'peachpuff3',
    'cornflowerblue')) # Plot the values in "Data", differentiating between the control and other models with color.
  text(Bar, 0, labels = c(names(Data)[1], Varnames_finder(names(Data)[2:length(names(Data))])), pos = 2, offset = 0,
    xpd = TRUE, srt = 90, cex = 3) # Format X tick labels.
  mtext('R-Squared', side = 2, line = 6, cex = 4, font = 2) } # Y axis label.

if (exists('Sites_modcompared')) { # If the variable "Sites_modcompared" exists...
  sapply(1:length(Sites_modcompared), ML_compareplotter) } # Apply the above function to each element in the variable.
