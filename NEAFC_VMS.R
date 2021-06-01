######################
##  NEAFC VMS DATA  ##
##    PROCESSING,   ##
##     May 2021     ##
#################nc###

rm(list=ls())

## libraries

library(ggplot2)
library(dplyr)
library(raster)
library(rasterize)
library(ncdf4)
library(data.table)
library(geosphere)


## data

work.year <-2020
#set working year to 2020

       vms <- read.csv(paste("data/", work.year, "/POSITIONS.csv", sep=""), header=T)
catch.data <- read.csv(paste("data/", work.year, "/CATCHES.csv",   sep=""), header=T) 
   vessels <- read.csv(paste("data/", work.year, "/VESSELS.csv",   sep=""), header=T)
## read in file of positions, catches, and vessel details
  
catch.data <- catch.data[catch.data$TM == "CAT",]   
## remove any catch on entry, exit or other messages
    
## in 2021, headers match, so no need to to change
#colnames(catch.data)[1] <- "CALLSIGN"
#vms$CALLSIGN <- paste(vms$CALLSIGN)
#vessels$CALLSIGN <- paste(vessels$CALLSIGN)
## make sure these are both factors

       vms$year <- as.numeric(substr(vms$DA, 1, 4))
       vms$time <- sprintf("%04d", vms$TI)
 vms$posix.time <- as.POSIXct(paste(substr(vms$DA, 1, 4), "-", substr(vms$DA, 5, 6), "-", substr(vms$DA, 7,8),
                                   " ", substr(vms$time, 1, 2), ":", substr(vms$time, 3,4), sep=""), format = "%Y-%m-%d %H:%M")
catch.data$time <-  sprintf("%04d", catch.data$TI)
catch.data$posix.time <- as.POSIXct(paste(substr(catch.data$DA, 1, 4), "-", substr(catch.data$DA, 5, 6), "-", substr(catch.data$DA, 7,8),
                                          " ", substr(catch.data$time, 1, 2), ":", substr(catch.data$time, 3,4), sep=""), 
                                    format = "%Y-%m-%d %H:%M")
catch.data$year <- substr(catch.data$DA, 1, 4)
## create a year field for catch and VMS data

##  1 VMS record outside 2020
##  271 catch records
## VMS record is from 23.55 31.12.2019, and matches a catch record (PRA in Barents)
## otherwise, dates are at various times in 2019 

write.csv(catch.data[catch.data$year != work.year,], "qc/catch_wrong_year.csv")
write.csv(vms[vms$year != work.year,], "qc/vms_wrong_year.csv")

vms <- vms[vms$year == work.year,]
catch.data <- catch.data[catch.data$year == work.year,]
## drop data outside 2020

    stored.vms <- vms
  stored.catch <- catch.data
stored.vessels <- vessels
## park data for now

dup.vms <- vms[duplicated(vms)==T,]
dup.catch.data <- catch.data[duplicated(catch.data)==T,]
dup.vessels <- vessels[duplicated(vessels)==T,]
## catch any duplicated records

write.csv(dup.vms, "qc/duplicated_vms.csv")
write.csv(dup.catch.data, "qc/duplicated_catch.csv")


vms <- vms[duplicated(vms)==F,]
catch.data <- catch.data[duplicated(catch.data)==F,]
vessels <- vessels[duplicated(vessels)==F,]
## remove any duplicate records


vms <- vms[order(vms$posix.time),]
## sort the data choronologically

unlinked.catch <- catch.data[catch.data$CALLSIGN %in% vms$CALLSIGN == FALSE,]
## 163 catch records without a corresponding vessel

unrecorded.vessels <- vessels[vessels$CALLSIGN %in% vms$CALLSIGN == FALSE,]
## 41 vessels not active in the VMS data


catch.data <- catch.data[catch.data$CALLSIGN %in% vms$CALLSIGN,]

vessels.not.fishing <- vessels[vessels$CALLSIGN %in% vms$CALLSIGN == FALSE,]
vessels.fishing <- vessels[vessels$CALLSIGN %in% vms$CALLSIGN==TRUE,]
unlinked.vms <- vms[vms$CALLSIGN %in% vessels$CALLSIGN ==FALSE,]

vessels <- vessels[vessels$CALLSIGN %in% vms$CALLSIGN,]
## drop vessels and catches from earlier than 2019

vms$Gear <- vessels$FISHING_GEAR[match(vms$CALLSIGN, vessels$CALLSIGN)]
## assign gear type to vessels using vessel registry table

vms$Gear <- paste(vms$Gear)
vms$Gear[vms$Gear == ""] <- "UNK"
## need to ask NEAFC if there is a difference between specified "NIL" and blank gear type


sort(tapply(rep(1, dim(vms)[1]), vms$Gear, sum))
## we still have 27k pings with no gear - will come back to these


vms.all <- vms
## park this for now too in case we screw up


callsigns <- unique(vms$CALLSIGN)
# 216 vessels fishing in 2019

sum(tapply(rep(1, dim(vms)[1]), vms$CALLSIGN, sum, na.rm=T)==0)
## no trips have zero pings this year, but in previous years some have one, or a just ENT and EXT value

sum(tapply(rep(1, dim(vms)[1]), vms$CALLSIGN, sum, na.rm=T)==1)
## that doesn't seem to be the case in 2020


vms <-vms[vms$TM %in% c("POS", "MAN"),]
## removes ENT & EXT records

## single.pings <- names(tapply(rep(1, dim(vms)[1]), vms$CALLSIGN, sum, na.rm=T)==1)[(tapply(rep(1, dim(vms)[1]), vms$CALLSIGN, sum, na.rm=T)==1)==TRUE]
## vms <- vms[vms$RID %in% single.pings == FALSE,]
## removes trips with a single ping (?!)
## don't have any of these in 2020

callsigns <- unique(vms$CALLSIGN)

sum(tapply(rep(1, dim(vms)[1]), vms$CALLSIGN, sum, na.rm=T)==0)

vms$CALLSIGN <- as.factor(vms$CALLSIGN)
vms <- data.frame(vms)

## this is a loop which orders VMS records chronologically by vessel, calculates the distance covered and speed

vms.out <- NULL

for (i in (1:length(callsigns))){
  
  
 temp.vms <- filter(vms, vms$CALLSIGN == callsigns[i])
  
 temp.vms <- temp.vms %>% 
    mutate(estimatedTimeDifHOurs = c(0, difftime(tail(posix.time, -1), head(posix.time, -1), units = "hours")))

    scratch <- distGeo(cbind(temp.vms$LG, temp.vms$LT))

    temp.vms$dist_nm <- c(0, scratch[-length(scratch)]/540)

    vms.out <- rbind(vms.out, temp.vms)
  
}

vms.out$ICES_Speed <- vms.out$SP / 10


hist(round(vms.out$SP/10,0), breaks = c(0:26))


vms.out <- vms.out[vms.out$estimatedTimeDifHOurs!=0,]
vms.out$estimatedTimeDifHOurs[is.na(vms.out$estimatedTimeDifHOurs)] <- mean(vms.out$estimatedTimeDifHOurs[!is.na(vms.out$estimatedTimeDifHOurs)])

vms.out$Actual.TDiff <- vms.out$estimatedTimeDifHOurs

vms.out$Actual.TDiff[vms.out$Actual.TDiff>2] <- 2


vms.all <- vms.out
## park this for now too

## Specify gear groupings
 trawl.gears <- c("OTB", "PTB", "TBS", "OTT")
static.gears <- c("LL", "LLS", "LLD", "GND", "GNS", "LNB")
    no.gears <- c("NULL", "NIL", "UNK")

vms.trawl  <- vms.all[vms.all$Gear %in% trawl.gears,]
vms.static <- vms.all[vms.all$Gear %in% static.gears,]
vms.nogear <- vms.all[vms.all$Gear %in% no.gears,]
## subset vms by gear type

dev.off()
par(mfrow=c(1,3))
hist(vms.trawl$ICES_Speed, breaks = 50, xlim=c(0,15), col="lightskyblue", main ="Bottom Trawls", xlab = "Speed (knots)")
hist(vms.static$ICES_Speed, breaks = 15, xlim=c(0,15), col="lightskyblue", main ="Static Gear", xlab = "Speed (knots)")
hist(vms.nogear$ICES_Speed, breaks = 15, xlim=c(0,15), col="lightskyblue", main ="No recorded gear", xlab = "Speed (knots)", )


## first we output vessels with no recorded gear

no.gear.points <- SpatialPointsDataFrame(data.frame(x = vms.nogear$LG[vms.nogear$ICES_Speed<=5], y = vms.nogear$LT[vms.nogear$ICES_Speed<=5]), data.frame(hours.fished = vms.nogear$Actual.TDiff[vms.nogear$ICES_Speed<=5]))
no.gear.raster <- raster(ncol = 700, nrow = 560, xmn = -42, ymn = 36, xmx = -7, ymx = 64)
no.gear.effort <- rasterize(no.gear.points, no.gear.raster, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(no.gear.effort, "results/no_gear_2021_hours_fished.asc", format = "CDF", overwrite=T)
## and writes it out in a format to be read in GIS software

## then we can do static gears

static.gear.points <- SpatialPointsDataFrame(data.frame(x = vms.static$LG[vms.static$ICES_Speed<=4], y = vms.static$LT[vms.static$ICES_Speed<=4]), data.frame(hours.fished = vms.static$Actual.TDiff[vms.static$ICES_Speed<=4]))
static.gear.raster <- raster(ncol = 700, nrow = 560, xmn = -42, ymn = 36, xmx = -7, ymx = 64)
static.gear.effort <- rasterize(static.gear.points, static.gear.raster, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(static.gear.effort, "results/static_gear_2021_hours_fished.asc", format = "CDF", overwrite = T)
## and writes it out in a format to be read in GIS software

## and finally the mobile gears

trawl.gear.points <- SpatialPointsDataFrame(data.frame(x = vms.trawl$LG[vms.trawl$ICES_Speed<=5], y = vms.trawl$LT[vms.trawl$ICES_Speed<=5]), data.frame(hours.fished = vms.trawl$Actual.TDiff[vms.trawl$ICES_Speed<=5]))
trawl.gear.raster <- raster(ncol = 700, nrow = 560, xmn = -42, ymn = 36, xmx = -7, ymx = 64)
trawl.gear.effort <- rasterize(trawl.gear.points, trawl.gear.raster, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(trawl.gear.effort, "results/trawl_gear_2021_hours_fished.asc", format = "CDF", overwrite = T)
## and writes it out in a format to be read in GIS software

### Now we do the same for NEAFC Regulatory Area 2
no.gear.raster.2 <- raster(ncol = 360, nrow = 300, xmn = -6, ymn = 62, xmx = 12, ymx = 77)
no.gear.effort.2 <- rasterize(no.gear.points, no.gear.raster.2, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(no.gear.effort.2, "results/no_gear_2021_hours_fished_RA2.asc", format = "CDF", overwrite=T)
## and writes it out in a format to be read in GIS software

trawl.gear.raster.2 <- raster(ncol = 360, nrow = 300, xmn = -6, ymn = 62, xmx = 12, ymx = 77)
trawl.gear.effort.2 <- rasterize(trawl.gear.points, trawl.gear.raster.2, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(trawl.gear.effort.2, "results/trawl_gear_2021_hours_fished_RA2.asc", format = "CDF", overwrite = T)
## and writes it out in a format to be read in GIS software


### ...and for NEAFC Regulatory Area 3. 
no.gear.raster.3 <- raster(ncol = 320, nrow = 160, xmn = 30, ymn = 70, xmx = 46, ymx = 78)
no.gear.effort.3 <- rasterize(no.gear.points, no.gear.raster.3, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(no.gear.effort.2, "results/no_gear_2021_hours_fished_RA3.asc", format = "CDF", overwrite=T)
## and writes it out in a format to be read in GIS software

trawl.gear.raster.3 <- raster(ncol = 320, nrow = 160, xmn = 30, ymn = 70, xmx = 46, ymx = 78)
trawl.gear.effort.3 <- rasterize(trawl.gear.points, trawl.gear.raster.3, field = "hours.fished", fun= sum)
## sums the hours fished data over the raster
writeRaster(trawl.gear.effort.3, "results/trawl_gear_2021_hours_fished_RA3.asc", format = "CDF", overwrite = T)
## and writes it out in a format to be read in GIS software



## The next piece of code groups consecutive pings at fishing speeds for mobile gears into putative "hauls"
## if you are not doing this straight after the above section, make sure vms.out <- vms.trawl

callsigns <- unique(paste(vms.trawl$CALLSIGN))

haul.id <- 1

haul.vms <- NULL
## set starting values and structures

for(i in (1:length(callsigns))){
  
  temp.vms <- vms.trawl[vms.trawl$CALLSIGN==callsigns[i],]
  temp.vms <- temp.vms[!is.na(temp.vms$ICES_Speed),]
  start.state  <- ifelse(temp.vms$ICES_Speed[1]<5, "fishing", "steaming")
  # determine if first ping is at fishing or steaming speed
  haul.id <- ifelse(start.state == "fishing", haul.id+1, haul.id)
  # if the vessel is fishing, add 1 to "HAUL_ID", otherwise, leave it along
temp.vms$Haul.No <- temp.vms$distance.covered <- NA
  
    temp.vms$Haul.No[1] <- ifelse(start.state=="fishing", haul.id, NA)
 
  for(j in (2:dim(temp.vms)[1])){
    # for each subsequent ping

    if(temp.vms$dist_nm[j]>10){
      haul.id <- haul.id + 1
      next 
    }  ## if a vessel has moved more than 10nm between pings, assume a new haul

    if(temp.vms$Actual.TDiff[j]>2){
      haul.id <- haul.id+1
      next
    } ## if the time between pings is greater than two hours, assume a new haul

    if(temp.vms$ICES_Speed[j] >5){temp.vms$Haul.No[j] <- NA}
      # if the vessel remains above 5 knots, no change to haul ID
      
    if(temp.vms$ICES_Speed[j] <= 1){
      temp.vms$Haul.No[j] <- NA
      haul.id <- haul.id + 1
    next
    }
      ## if the vessel is not moving, ignore
    

    if(temp.vms$ICES_Speed[j] <=5 & temp.vms$ICES_Speed[j] >= 1 & temp.vms$ICES_Speed[j-1] > 5){
      haul.id <- haul.id + 1
      temp.vms$Haul.No[j] <- haul.id
    next
    }     # if the vessel slows down to between 1 and 5 knots, assume a new haul and assign it to the ping


    if(temp.vms$ICES_Speed[j] <= 5 & temp.vms$ICES_Speed[j] >=1 & 
       temp.vms$ICES_Speed[j-1] <= 5 & temp.vms$ICES_Speed[j-1] >=1) {
      temp.vms$Haul.No[j] <- haul.id
      next
    }
    # if the vessel continues to move at 1-5 knots, it remains fishing and no change to haul ID

    if(temp.vms$ICES_Speed[j] > 5 & temp.vms$ICES_Speed[j-1] <= 5 &
       temp.vms$ICES_Speed[j-1] >=1 ) {
      haul.id <- haul.id + 1
      temp.vms$Haul.No[j] <- NA
      next
      }

  } 
  haul.vms <- rbind(haul.vms, temp.vms)
  # binds temporary data to output
}


hauls.out <- haul.vms[!is.na(haul.vms$Haul.No),]
## drop rows where there isn't fishing activity

for(i in (1:length(unique(hauls.out$Haul.No)))){
  hauls.out$Sequence[hauls.out$Haul.No == unique(hauls.out$Haul.No)[i]] <- c(1:sum(hauls.out$Haul.No==unique(hauls.out$Haul.No)[i]))
}



hauls.out

hauls.out <-  hauls.out[hauls.out$Haul.No %in% unique(hauls.out$Haul.No)[tapply(rep(1, dim(hauls.out)[1]), hauls.out$Haul.No, sum)>1],]
## drop hauls where we only have one point - no way to plot lines

hauls.out <- hauls.out[,colnames(hauls.out) %in% c("LT", "LG", "Haul.No", "Sequence", "distance.covered")]
## drops other columns for outputting


write.csv(hauls.out, "C:/Work/WGSFD/2021/OTB_TOWS_2021.csv", row.names=TRUE)
