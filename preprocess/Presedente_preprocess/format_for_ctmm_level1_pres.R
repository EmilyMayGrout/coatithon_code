#occurrence distribution plot for Presidente group

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level1.RData"
id_file <- 'presedente_coati_ids.RData' 

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

Ardern <- data.frame(cbind(xs[1,],ys[1,]))
Ardern$ts <- as.POSIXct(ts)

Torrijos <- data.frame(cbind(xs[19,],ys[19,]))
Torrijos$ts <- as.POSIXct(ts)

Wildflower <- data.frame(cbind(xs[21,],ys[21,]))
Wildflower$ts <- as.POSIXct(ts)

#........................
library(ctmm)
library(tidyverse)
library(sf)
library(ggmap)

ard <- Ardern %>% 
  mutate(individual.local.identifier = "Ardern") %>% 
  rename(X = X1,
         Y = X2,
         timestamp = ts) %>% 
  drop_na()

torr <- Torrijos %>% 
  mutate(individual.local.identifier = "Torrijos") %>% 
  rename(X = X1,
         Y = X2,
         timestamp = ts) %>% 
  drop_na()

wild <- Wildflower %>% 
  mutate(individual.local.identifier = "Wildflower") %>% 
  rename(X = X1,
         Y = X2,
         timestamp = ts) %>% 
  drop_na()


DATA <- rbind(ard, torr, wild) %>% 
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
plot(DATA, col = c("aquamarine3"), main = "Location Data")  

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
par(mfrow = c(2,2))
# plot empirical variogram with best model
plot(SVF[1], CTMM = FIT[1], main = "Ardern")
plot(SVF[2], CTMM = FIT[2], main = "Torrijos")
#plot(SVF[3], CTMM = FIT[3], main = "Wildflower")

# summary
summary(FIT[[1]])
summary(FIT[[2]])
#summary(FIT[[3]])
summary(UD[[1]]) #for the area
summary(UD[[2]])
#summary(UD[[3]])
# plot home ranges
plot(DATA[1],
     UD=UD[1],
     col = "#e76f51",
     col.DF="#f4a261",
     col.grid = NA,
     main = "Ardern")
plot(DATA[2],
     UD=UD[2],
     col = "#264653",
     col.DF="#2a9d8f",
     col.grid = NA,
     main = "Torrijos")
plot(DATA[3],
     UD=UD[3],
     col = "#e76f51",
     col.DF="#f4a261",
     col.grid = NA,
     main = "Wildflower")

saveRDS(UD, "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/AKDEs_presedente.rds")

#to get the area of one individual
area_ard <- summary(UD[[1]], units = FALSE)$CI[2]/1000000
area_torr <- summary(UD[[2]], units = FALSE)$CI[2]/1000000
#area_wild <- summary(UD[[3]], units = FALSE)$CI[2]/1000000

#looking at one individual on an interactive map
UD_sf <- as.sf(UD[[1]])
mapview(UD_sf)

#get the area and confidence intervals
#if only want to make a plot for group members, need to remove wildflower from the DATA in line 55

UD_sf <- UD %>%
  purrr::map(ctmm::as.sf) %>% 
  reduce(rbind) %>%
  filter(grepl("est", name)) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
  mutate(coati = gsub( " .*$", "", name))

CI_sf <- UD %>%
  purrr::map(ctmm::as.sf) %>% 
  reduce(rbind) %>%
  filter(!grepl("est", name)) %>% 
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
  mutate(coati = gsub( " .*$", "", name))

DATA_sf <- DATA %>%
  purrr::map(ctmm::as.sf) %>% 
  reduce(rbind) %>%
  st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
  mutate(coati = gsub( " .*$", "", identity))

#need to register a google key for the map to work
map = get_map(location = c(lon =-79.837603, lat=9.168051), zoom=16, maptype="satellite")

#if we only want to look at one individual, need to subset data: e.g UD_sf[UD_sf$name[3],] but for each one

gg <- ggmap(map) + 
  geom_sf(data = UD_sf, 
          #color = "#d3436e",
          aes(fill = coati),
          inherit.aes = FALSE, 
          alpha = 0.3, 
          size = 15,
          color = "white") +
  geom_sf(data = CI_sf, 
          fill = NA,
          inherit.aes = FALSE, 
          alpha = 0.2, 
          size = 10,
          color = "white",
          lty ="dotted") +
  geom_sf(data = DATA_sf,
          color = "lavenderblush",
          shape = 1,
          size = 0.1,
          aes(color = coati),
          inherit.aes = FALSE) +
  scale_fill_viridis_d(end = 0.5) +
  #facet_wrap(~coati, ncol = 2)+
# theme(axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())
NULL

gg

#ggsave(filename = paste0(plot_dir, 'overlap', '.png'), plot = gg, width = 8, height = 8, dpi = 500)


#calculate distance travelled between subgroups - using distance between Estrella and Planeta


dist_between_inds  <- cbind(Ardern, Torrijos)
colnames(dist_between_inds) <- c("ard_eastings", "ard_northings", "ard_ts", "tor_eastings", "tor_northings", "tor_ts")
dist_between_inds$x_dist <- (dist_between_inds$ard_eastings - dist_between_inds$tor_eastings)
dist_between_inds$y_dist <- (dist_between_inds$ard_northings - dist_between_inds$tor_northings)
dist_between_inds$dist <- sqrt(((dist_between_inds$x_dist)*(dist_between_inds$x_dist))+((dist_between_inds$y_dist)*(dist_between_inds$y_dist)))


#png(height = 800, width = 1000, units = 'px', filename = paste0(plot_dir,"dist_between_ard_tor.png"))

par(mar=c(5.1,5.1,4.1,2.1))# bottom, left, top and right 
hist(dist_between_inds$dist, breaks = 200, col = "lightblue", main = " ", xlab = "Distance between Ardern and Torrijos (m)", cex.lab = 2, cex.axis= 1.5)

dev.off()

#-------------------------------------------------------------------------------------

#calculate the daily travel distance for Ardern

# Create an sf data frame for one ind
sf_data <- DATA$Wildflower

dpl <- sf_data %>%
  as.sf() %>% 
  mutate(date = as_date(timestamp)) %>% 
  group_by(date) %>% 
  mutate(lead = geometry[row_number() + 1],
         dist = st_distance(geometry, lead, by_element = T))

dpl_sum <- dpl %>% 
  mutate(dist = as.numeric(dist)) %>%
  na.omit() %>% 
  group_by(date) %>% 
  dplyr::summarise(dpl = sum(dist))











