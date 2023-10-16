#occurance distribution plot for Galaxy group

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/'
gps_file <- "galaxy_xy_10min_level1.RData"
id_file <- 'coati_ids.RData' 

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

Lucero <- data.frame(cbind(xs[4,],ys[4,]))
Lucero$ts <- as.POSIXct(ts)

Estrella <- data.frame(cbind(xs[2,],ys[2,]))
Estrella$ts <- as.POSIXct(ts)


#........................
library(ctmm)
library(tidyverse)
library(sf)

luc <- Lucero %>% 
  mutate(individual.local.identifier = "Lucero") %>% 
  rename(X = X1,
         Y = X2,
         timestamp = ts) %>% 
  drop_na()

est <- Estrella %>% 
  mutate(individual.local.identifier = "Estrella") %>% 
  rename(X = X1,
         Y = X2,
         timestamp = ts) %>% 
  drop_na()

DATA <- rbind(luc,est) %>% 
  st_as_sf(crs="+proj=utm +zone=17 +north +datum=WGS84 +units=m", 
           coords = c("X", "Y")) %>% 
  st_transform(crs = st_crs("+proj=longlat +datum=WGS84")) %>% 
  mutate(location.long = st_coordinates(geometry)[, 1],
         location.lat = st_coordinates(geometry)[, 2]) %>% 
  st_drop_geometry() %>% 
  dplyr::select(individual.local.identifier,
                timestamp,
                location.long,
                location.lat) %>% 
  as.telemetry(projection = "+proj=utm +zone=17 +north +datum=WGS84 +units=m +no_defs
↪ +ellps=WGS84 +towgs84=0,0,0")

# plot location data
plot(DATA, col = c("aquamarine3", "orange2" ), main = "Location Data")  

# variogram
UD <- FIT <- SVF <- list()
for(i in 1:length(DATA)){
  SVF[[i]] <- variogram(DATA[[i]])
  GUESS <- ctmm.guess(DATA[[i]],
                      interactive=FALSE,
                      variogram = SVF[[i]])
  FIT[[i]] <- ctmm.select(DATA[[i]],
                           GUESS,
                           trace=2)
  UD[[i]] <- akde(DATA[[i]],
                   FIT[[i]],
                   grid=list(dr=10, align.to.origin=TRUE))
}

names(SVF) <- names(FIT) <- names(UD) <- names(DATA)

# plot variograms and ctmms
par(mfrow = c(1,2))
# plot empirical variogram with best model
plot(SVF[1], CTMM = FIT[1], main = "Lucero")
plot(SVF[2], CTMM = FIT[2], main = "Estrella")

# summary
summary(FIT[[1]])
summary(FIT[[2]])
summary(UD[[1]])
summary(UD[[2]])

# plot home ranges
plot(DATA[1],
     UD=UD[1],
     col = "#e76f51",
     col.DF="#f4a261",
     col.grid = NA,
     main = "Lucero")
plot(DATA[1],
     UD=UD[1],
     col = "#264653",
     col.DF="#2a9d8f",
     col.grid = NA,
     main = "Estrella")

saveRDS(UD, "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/AKDEs_galaxy.rds")
