rm(list=ls())
library(icesVMS)
library(icesVocab)
library(sf)
library(data.table)
library(dplyr)
library(tidyr)
library(vmstools)
sf_use_s2(FALSE)

## Process a polygon covering the coastline of Div. 3a
data(ICESareas)
ICESareas <- ICESareas %>%
  dplyr::filter(Area_27 %in% c("3.a.20", "3.a.21"))
ICESareas <- st_union(ICESareas)

## Download the c-squares in 3a from ICES
kat.csq <- get_csquare(ices_area = "3a", convert2sf = TRUE)

# Ensure both objects are in WGS84
st_crs(kat.csq) <- st_crs(4326)
st_crs(ICESareas) <- st_crs(4326)

# Crop kat.csq to ICESareas
kat.csq_cropped <- st_intersection(kat.csq, ICESareas)

# Function to get an appropriate UTM zone based on longitude
get_utm_zone <- function(lon) {
  zone <- floor((lon + 180) / 6) + 1
  return(paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs"))
}

# Calculate centroid longitude for the entire dataset
centroid_lon <- st_coordinates(st_centroid(st_as_sfc(st_bbox(kat.csq_cropped))))[1]

# Get appropriate UTM CRS
utm_crs <- get_utm_zone(centroid_lon)

# Project to UTM, calculate area, and project back to WGS84
kat.csq_result <- kat.csq_cropped %>%
  st_transform(utm_crs) %>%
  mutate(area_km2 = as.numeric(st_area(.) / 1e6)) %>%  # Convert m^2 to km^2
  st_transform(4326)

# Keep only the original fields and the new area
kat.csq_result <- kat.csq_result %>%
  select(names(kat.csq), area_km2)

# Load the specific layer from the geodatabase
gdb_path <- "C:/Work/WGSFD/ICES-VMS-and-Logbook-Data-Call_refactor/data/eusm.rds"
eusm <- readRDS(gdb_path)

# Keep only the "MSFD_BBHT" variable and join with codelist
codelist <- getCodeList("BroadHabitat")
eusm <- eusm %>% 
  dplyr::select(MSFD_BBHT) %>%
  left_join(codelist, by = c("MSFD_BBHT" = "Description")) %>%
  dplyr::mutate(MSFD_BBHT = dplyr::case_when(
    MSFD_BBHT == "Na" ~ "",
    is.na(Key) ~ "",
    TRUE ~ Key
  )) %>%
  dplyr::select(MSFD_BBHT)

## set up a bounding box slightly larger than the c-square object
bbox <- st_bbox(c(xmin = min(kat.csq$lon)-0.1, 
                        ymin = min(kat.csq$lat)-0.1, 
                        xmax = max(kat.csq$lon + 0.1), 
                        ymax = max(kat.csq$lat + 0.1)), crs = 4326)

bbox_wgs84 <- st_transform(st_as_sfc(bbox), st_crs(eusm))

## Crop the shapefile to the bounding box
eusm_filtered <- st_crop(eusm, bbox_wgs84)

rm(eusm)

# Create a grid of points covering Div. 3a, 100 per c-square
grid_size <- 0.05
points_grid <- st_make_grid(ICESareas, cellsize = c(grid_size/3, grid_size/3))

points_df <- points_grid %>%
  st_cast("POINT") %>%
  st_cast("MULTIPOINT") %>%
  st_cast("POINT") %>%
  st_as_sf()

# Convert the grid to an sf object
points_sf <- st_centroid(points_grid)

points_sf <- st_as_sf(data.frame(geometry = points_sf))

# Add coordinates as columns if needed
points_sf <- points_sf %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2])

# Add a unique identifier if needed
points_sf <- points_sf %>%
  mutate(id = row_number())

# Now perform the spatial join
points_habitat <- st_join(points_sf, eusm_filtered)

points_habitat$CSquare <- vmstools::CSquare(points_habitat$lon, points_habitat$lat, 0.05)

st_crs(points_habitat) <- st_crs(kat.csq_result)

# Clean MSFD_BBHT column: replace empty strings, NA, and whitespace-only values with "Unknown"
points_habitat <- points_habitat %>%
  mutate(MSFD_BBHT = ifelse(is.na(MSFD_BBHT) | MSFD_BBHT == "" | trimws(MSFD_BBHT) == "", "Unknown", MSFD_BBHT))

# Count points of each habitat type within each CSquare
habitat_counts <- points_habitat %>%
  st_drop_geometry() %>%  # Remove geometry to speed up grouping
  group_by(CSquare) %>%
  mutate(total_points = n()) %>%
  group_by(CSquare, MSFD_BBHT, total_points) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total_points) %>%
  select(-count, -total_points)

# Check for any remaining problematic values
print(unique(habitat_counts$MSFD_BBHT))

# If there are still empty strings, replace them with "Unknown"
habitat_counts <- habitat_counts %>%
  mutate(MSFD_BBHT = ifelse(MSFD_BBHT == "", "Unknown", MSFD_BBHT))

# Pivot the data to create columns for each habitat type
habitat_proportions_wide <- habitat_counts %>%
  pivot_wider(
    id_cols = CSquare,
    names_from = MSFD_BBHT,
    values_from = proportion,
    values_fill = list(proportion = 0)
  )

# Join the proportions back to the original polygon data
kat.csq_result_with_areas <- kat.csq_result %>%
  left_join(habitat_proportions_wide, by = c("c_square" = "CSquare"))

# Replace NA values with 0 for the new habitat columns
habitat_columns <- setdiff(names(habitat_proportions_wide), "CSquare")
kat.csq_result_with_areas[habitat_columns] <- lapply(kat.csq_result_with_areas[habitat_columns], function(x) replace(x, is.na(x), 0))

# Multiply proportions by area_km2 to get estimated habitat areas
for (col in habitat_columns) {
  kat.csq_result_with_areas[[col]] <- kat.csq_result_with_areas[[col]] * kat.csq_result_with_areas$area_km2
}

# Ensure the result is still an sf object (it should be, but just to be safe)
kat.csq_result_with_areas <- st_as_sf(kat.csq_result_with_areas)

# Preview the result
print(head(kat.csq_result_with_areas))

                                                     
