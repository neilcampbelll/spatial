#-------------------------------------------------------------------------------
#
# Script to extract and process VMS and logbook data for ICES VMS data call
#
# By: Niels Hintzen, Katell Hamon, Marcel Machiels#
# Code by: Niels Hintzen
# Contact: niels.hintzen@wur.nl
#
# Date: 25-Jan-2017
# Update Date: 29-Jan-2019 ; Updated by: Roi Martinez
# Update Date: 04-Feb-2020 ; Updated by: Colin Millar
# Update Date: 07-Feb 2020 ; Updated by: Neil Campbell
# Client: ICES
#-------------------------------------------------------------------------------

#--------------------READ ME----------------------------------------------------
# The following script is a proposed workflow example to processes the ICES
# VMS datacall request. It is not an exact template to be applied to data from
# every member state and needs to be adjusted according to the data availability
# and needs of every member state.
#-------------------------------------------------------------------------------


#- Clear workspace
rm(list=ls())

library(vmstools) #- download from www.vmstools.org
library(Matrix)   #- available on CRAN
library(ggplot2)  #- available on CRAN
library(dplyr)    #- available on CRAN
library(sp)
library(doBy)
library(mixtools)
library(tidyr)
library(glue)
library(gt)
library(raster)
library(sf)
library(data.table)
 
#- Settings paths
codePath  <- "R"          #Location where you store R scripts
dataPath  <- "Data"       #Location where you store tacsat (VMS) and eflalo (logbook) data
outPath   <- "Results"    #Location where you want to store the results

#- Setting specific thresholds
spThres       <- 20   #Maximum speed threshold in analyses in nm
intThres      <- 5    #Minimum difference in time interval in minutes to prevent pseudo duplicates
intvThres     <- 240  #Maximum difference in time interval in minutes to prevent intervals being too large to be realistic
lanThres      <- 1.5  #Maximum difference in log10-transformed sorted weights

#- Re-run all years as we have new field for no. vessels
yearsToSubmit <- c(2018, 2022)

#- Set the gear names for which automatic fishing activity is wanted
#  It is important to fill out the gears you want to apply auto detection for
autoDetectionGears        <- c("DRB_MOL", "FPO_CRU", "FPO_DEF", "FPO_MOL", "GNS_DEF", "GNS_SPF",
                               "GTR_DEF", "LLD_LPF", "LLS_DEF", "LLS_DWS", "MIS_MIS",
                                "OTB_CRU", "OTB_DEF", "OTB_MCD", "OTB_MCF")

#- Decide if you want to visualy analyse speed-histograms to identify fishing activity
#  peaks or have prior knowledge and use the template provided around lines 380 below
visualInspection          <- TRUE

#- Specify how landings should be distributed over the VMS pings: By day, ICES rectangle, trip basis or otherwise
linkEflaloTacsat          <- c("day","ICESrectangle","trip")
# other options
# linkEflaloTacsat          <- c("day","ICESrectangle","trip")
# linkEflaloTacsat          <- c("ICESrectangle","trip")
# linkEflaloTacsat          <- c("day","trip")
# linkEflaloTacsat          <- c("trip")

