library(dplyr)
library(icesVMS)
library(vmstools)
library(sf)

for(i in 2009:2023){  
  
  vms <- readRDS(paste0("data/VMS", i, "_with_swept_area.RDS"))
  vms$anonVessels[vms$uniqueVessels >= 3] <- "not_required"
  ## in "normal times" this line would be replaced with a call to the VMS database
  
  # Filter data 
  vms_filtered <- vms %>%
    dplyr::select(country, year, month, cSquare, leMetLevel6, fishingHours, 
                  kwFishinghours, totweight, totvalue, uniqueVessels, anonVessels) %>%
    filter(country %in% c("BE", "DE", "NL", "DK"), 
           leMetLevel6 %in% c("TBB_CRU_16-31_0_0", "TBB_DEF_16-31_0_0"))
  
  # Function to process and save data
  process_and_save <- function(data, grouping_vars, filename_suffix) {
    # Check if data is empty
    if(nrow(data) == 0) {
      cat("No data found for", filename_suffix, "in year", i, "- skipping\n")
      return(NULL)
    }
    
    temp <- data %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarise(
        fishingHours = round(sum(fishingHours), 3),
           totweight = round(sum(totweight), 3),
                LPUE = round(sum(totweight) / sum(fishingHours), 3),
            LPUE_KWH = round(sum(totweight) / sum(kwFishinghours), 3),
       
         # Check if any record already had "not_required" status
        has_not_required = any(anonVessels == "not_required"),
        
        # Collect all individual vessel IDs (excluding "not_required")
        all_vessel_ids = list(unique(unlist(strsplit(anonVessels[anonVessels != "not_required"], ";")))),
        .groups = 'drop'
      ) %>%
      mutate(
        # Count unique vessels
        vessel_count = ifelse(
          has_not_required,
          3, # If any record had "not_required", assume 3+ vessels total
          lengths(all_vessel_ids)
        ),
      
          # Create mask - TRUE if needs masking (< 3 vessels), FALSE otherwise  
        Mask = vessel_count < 3,
        
        # Final vessel count and anonymisation status
        vessels = vessel_count,
        anonVessels = ifelse(
          vessel_count >= 3, 
          "not_required", 
          sapply(all_vessel_ids, paste, collapse = ";")
        )
      ) %>%
      select(-has_not_required, -all_vessel_ids, -vessel_count, -vessels, -anonVessels)
    
      
    # convert cSquares to geometry
      coords <- vmstools::CSquare2LonLat(temp$cSquare, 0.05)
      temp$geometry <- wkt_csquare(coords$SI_LATI, coords$SI_LON)
    
    # convert the data frame to a spatial object
    temp <- st_as_sf(temp, wkt = "geometry", crs = 4236)
    
    # save as text file
    temp_df <- st_drop_geometry(temp)
    temp_df$geometry <- st_as_text(temp$geometry)
    write.table(temp_df, paste0("output/text_files/CRAN_DATA_", filename_suffix, "_", i, ".wkt"), 
                sep = ";", row.names = FALSE)
  }
  
  # 1. TBB_DEF_16-31_0_0 only
  vms_def <- vms_filtered %>%
    filter(leMetLevel6 == "TBB_DEF_16-31_0_0")
  process_and_save(vms_def, c("cSquare"), "DEF_only")
  
  # 2. TBB_CRU_16-31_0_0 only  
  vms_cru <- vms_filtered %>%
    filter(leMetLevel6 == "TBB_CRU_16-31_0_0")
  process_and_save(vms_cru, c("cSquare"), "CRU_only")
  
  # 3. Both metiers combined
  process_and_save(vms_filtered, c("cSquare"), "combined")
}
