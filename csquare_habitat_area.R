rm(list=ls())

library(icesVMS)
library(sf)
library(dplyr)
library(units)
library(lwgeom)
library(vmstools) 

sf_use_s2(FALSE) # this is needed to stop the habitat layer intersecting with itself and causing problems

# Load the specific layer from the geodatabase
kat.csq <- get_csquare(ices_area = "3a", convert2sf = TRUE)

gdb_path <- "C:/Work/WGSFD/ICES-VMS-and-Logbook-Data-Call_refactor/Results/EUSeaMap_2023.gdb"
eusm <- readRDS(gdb_path)

# Keep only the "MSFD_BBHT" variable
eusm <- eusm %>% dplyr::select(MSFD_BBHT)
codelist <- getCodeList("BroadHabitat")

# Replace full names with abbreviations and handle non-matching habitat types
eusm <- eusm %>%
  left_join(codelist, by = c("MSFD_BBHT" = "Description")) %>%
  dplyr::mutate(MSFD_BBHT = dplyr::case_when(
    MSFD_BBHT == "Na" ~ "",
    is.na(Key) ~ "",
    TRUE ~ Key
  )) %>%
  dplyr::select(MSFD_BBHT)

# To speed the function up, create bounding box (0.05 x 0.05 degree)
# and filter the habitat layer to this first
bbox_wgs84 <- st_bbox(c(xmin = min(kat.csq$lon)-0.1, 
                    ymin = min(kat.csq$lat)-0.1, 
                    xmax = max(kat.csq$lon + 0.1), 
                    ymax = max(kat.csq$lat + 0.1)), crs = 4326)
  
                          
# Convert bbox to an sf polygon
bbox_poly_wgs84 <- st_as_sfc(bbox_wgs84)
                          
# Transform bbox polygon to the CRS of eusm (which is Web Mercator, EPSG:3857)
bbox_poly_mercator <- st_transform(bbox_poly_wgs84, st_crs(eusm))
                          
# Filter the data
eusm_filtered <- eusm[st_intersects(eusm, bbox_poly_mercator, sparse = FALSE)[,1], ]

# Create a function for processing each c-square
process_csquare <- function(wkt, csquare) {
  csq <- st_as_sf(wkt)
  
  # Calculate the area of this specific c-square
  csq_area_km2 <- as.numeric(units::set_units(st_area(csq), km^2))
  
  set.seed(123)  # for reproducibility
  points <- st_sample(csq, size = 1000, type = "regular")
  points_sf_transformed <- st_transform(st_set_crs(st_as_sf(points), 4326), st_crs(eusm_filtered))
  
  sampled_points <- st_join(points_sf_transformed, eusm_filtered)
  
  habitat_counts <- sampled_points %>%
    st_drop_geometry() %>%
    group_by(MSFD_BBHT) %>%
    summarise(count = n(), .groups = 'drop')
  
  results <- habitat_counts %>%
    mutate(
      area_estimate_km2 = (count / 1000) * csq_area_km2,
      csquare = csquare
    ) %>%
    dplyr::select(csquare, MSFD_BBHT, area_estimate_km2)
  
  return(results)
}

# Use lapply for processing
results_list <- lapply(1:nrow(kat.csq), function(i) {
  result <- process_csquare(kat.csq$wkt[i], kat.csq$c_square[i])
  if (i %% 100 == 0) print(paste("Processed", i, "of", nrow(kat.csq)))
  return(result)
})

# Combine results
results.out <- do.call(rbind, results_list)

# Convert to data.table for faster operations if needed
results.out <- as.data.table(results.out)

# Print the results
print(results.out)
