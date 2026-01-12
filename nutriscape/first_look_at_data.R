#visualising the BCI nutriscapes data - running the fission-fusion code to see the social associations

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 50

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level0/'
gps_file <- "presidente2025_xy_10min_level0.RData"
id_file <- 'presidente2025_coati_ids.RData'

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vioplot)
library(plotly)
library(rnaturalearth)
library(ggplot2)
library(gganimate)
library(ggmap)
library(ggspatial)
library(sp)
library(av)
library(cocomo)
library(patchwork)
library(sf)
library(ctmm)
library(reshape2)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)


#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds-1)
#low number because Ardern was on a lower res schedule, so minus 1 to ignore her

green <- rgb(120, 170, 80, maxColorValue = 255)

darkgreen <- rgb(120, 160, 80, maxColorValue = 255)

#-------------------------

#how much data are there per individual
xs_mat <- xs
rownames(xs_mat) <- coati_ids$name
long_xs <- melt(xs_mat, varnames = c("Row", "Col"), value.name = "value")

# Create a presence column: TRUE if value is present, FALSE if NA
long_xs$presence <- !is.na(long_xs$value)

ggplot(long_xs, aes(x = Col, y = Row)) +
  geom_point(aes(color = presence), size = 4, shape = 15) +
  scale_color_manual(values = c("grey90", "slateblue3"), labels = c("Absent", "Present")) +
  labs(x = "", y = "", color = "Data Present") +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

#------------------------

subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs], main = "", xlab =  'Number of subgroups (radius = 50 m)', col = green, border = "darkolivegreen3", breaks = length(unique(subgroup_data$n_subgroups[all_tracked_idxs])))

subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], main = "")
hist(subgroup_counts[,s3], main = "")

#-------------------------
png(height = 1400, width = 1200, units = 'px', res = 144, filename = paste0(plot_dir,'subgroups_diffradius.png'))

par(mfrow=c(3,2)) #(bottom, left, top, right)

for (i in 1:length(Rs)){
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = "Number of subgroups", col = "slateblue2")
}

dev.off()


#-------------------------

R <- 50

#Figure 3a: which individuals tend to be in the same subgroup
#subsample to times when all inds are tracked
xs_all <- xs[, all_tracked_idxs]
ys_all <- ys[, all_tracked_idxs]

subgroup_data <- get_subgroup_data(xs, ys, R)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
new_order <- c(1:16)
new_order <-c(1,4,7,15,6,13,16,9,2,10,5,11,12,14,3,8)
ffnet_reorder <- ff_net[new_order, new_order]
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'withingroup_network_difforder_', R, 'm.png'))

visualize_network_matrix_presedente(ffnet_reorder, coati_ids[new_order,])

dev.off()


png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'subgroupsize_', R, 'm.png'))

h <- hist(subgroup_data$subgroup_counts[, all_tracked_idxs],
          main = "", xlab = "Subgroup size", col = "slateblue3",
          xaxt = "n")  # suppress default x-axis

# Manually add x-axis with labels centered under the bars
axis(side = 1, at = h$mids, labels = 1:14)

dev.off()


#-----------------------------------------------------
#visualise movement trajectories - put function in coati_function_library_V1
pres_df <- matrix_to_df(xs, ys, ts)



#going through each day, extracting the mean lat and lon locations to zoom in the map and plot the movement trajectories 

#remove Ardern as not enough GPS points
coati_ids$ID <- c(1:16)
pres_df <- pres_df[!pres_df$ID == 1, ]
pres_df <- pres_df %>%
  left_join(coati_ids[, c("ID", "name")], by = "ID")

point_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928', "black", "white", "grey", "yellow", "darkblue")

run_code <- F

if(run_code == T){
for(i in 1:length(unique(as_date(pres_df$datetime)))){

date_i <- unique(as_date(pres_df$datetime))[i]
date_i_df <- pres_df[as_date(pres_df$datetime) == date_i,]
mean_lon <- mean(date_i_df$lon)
mean_lat <- mean(date_i_df$lat)

map = get_map(location = c(lon = mean_lon, lat = mean_lat), zoom=16, maptype="satellite")

g <- ggmap(map)+ 
  geom_path(data=date_i_df, aes(x = lon, y=lat, group= ID, colour= as.factor(name)))+
  geom_point(data=date_i_df, aes(x = lon, y=lat, group=ID, colour=as.factor(name)), size=4)+
  scale_color_manual(values = point_colors) +
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))

g <- g + transition_reveal(date_i_df$datetime, keep_last = FALSE)+
  ease_aes('linear')+
  labs(title = "Time: {format(as.POSIXct(frame_along, tz = 'UTC'), '%Y-%m-%d %H:%M:%S')}")+
  shadow_wake(wake_length = 0.05, alpha = FALSE)

#animate and save the video
vid <- animate(g, renderer = av_renderer(), height = 1000, width = 1000, fps = 10, duration = 40)

anim_save(filename = paste0(plot_dir, "animations/movement_", date_i, ".mp4"), animation = vid)
 }
}


#plot all raw data to see outliers
plot(pres_df$UTM_X, pres_df$UTM_Y, type = "p", asp = 1, cex = 0.1)


#clean the data to remove GPS outliers
breaks <- (which(diff(as_date(ts)) == 1)) + 1 #adding 1 as it gets the last each date and we need the index of the first time point for each day

#find which max_isolated_point_distance catches the errors but not the real data
mipd <- seq(from = 50, to = 1000, by = 10)
na_counts <- as.data.frame(mipd)
na_counts$na <- NA

if(run_code == T){
for(i in 1:length(mipd)){
clean <- preprocess_gps_level0_to_level1(xs = xs, ys = ys, timestamps = ts, ids = coati_ids, breaks = breaks, max_isolated_point_dist = mipd[i], max_dist_percentile = 0.98, verbose = F)
na_counts$na[i] <- sum(is.na(clean$xs))
}

plot(na_counts$mipd, na_counts$na, type = "l", xlab = "max_isolated_point_dist", ylab = "NA count")
}

clean <- preprocess_gps_level0_to_level1(xs = xs, ys = ys, timestamps = ts, ids = coati_ids, breaks = breaks, max_isolated_point_dist = 85, max_dist_percentile = 0.98, max_speed_percentile = 0.999, verbose = F)


#use function to convert matrix to dataframe for plotting the clean data
clean_UTM_df <- matrix_to_df(xs = clean$xs, ys = clean$ys, ts = clean$timestamps)

clean_UTM_df <- clean_UTM_df %>%
  left_join(coati_ids[, c("ID", "name")], by = "ID")

map = get_map(location = c(lon = mean(clean_UTM_df$lon), lat = mean(clean_UTM_df$lat)), zoom=16, maptype="satellite")

p <- ggmap(map)+ 
  geom_point(data=clean_UTM_df, aes(x = lon, y=lat),col="white", size=0.2, alpha = 0.15)+
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))+
  ggtitle("max_isolated_point_dist = 85m")
p


#ggsave(filename = paste0(plot_dir, 'potentialhotspots.png'), plot = p, width = 8, height = 8, dpi = 300)


#------------------------------------------------------

#look at the dipteryx tree crowns

dipx <- st_read("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")
dipx <- st_transform(dipx, crs = 4326)

# Get bounding box of shapefile
bbox <- st_bbox(dipx)

# Get base map
basemap <- get_map(
  location = c(lon = mean(clean_UTM_df$lon), 
               lat = mean(clean_UTM_df$lat)),
  zoom = 16,   # adjust zoom level as needed
  maptype = "satellite"  # or "satellite", "roadmap", etc.
)


# First, fix invalid geometries
dipx <- st_make_valid(dipx)

# Then calculate centroids
dipx_centroids <- st_centroid(dipx)

dipx_centroids <- dipx_centroids %>%
  st_transform(4326) %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  )

# Plot shapefile over basemap
hotspots <- ggmap(basemap) +
  theme_minimal()+ 
  geom_point(data=clean_UTM_df, aes(x = lon, y=lat),col="white", size=0.2, alpha = 0.15)+
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))+
  geom_sf(data = dipx, inherit.aes = FALSE, fill = NA, color = "magenta") +
  geom_text(data = dipx_centroids,
            aes(x = X, y = Y, label = newID),
            size = 1.5, color = "red")

#ggsave(filename = paste0(plot_dir, 'dipteryx_hotspots.png'), plot = hotspots, width = 10, height = 10, dpi = 300)


#looking at csv Kate made
crowns <- read.csv("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/CoatiTreeCrowns_2024-12-10_2025-02-16.csv")


#add crowns to animation with the cleaned GPS data
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level1/'

#remove Ardern from visualisations as she's on a low res
clean_UTM_df <- clean_UTM_df[clean_UTM_df$name != "Ardern",]

if(run_code == T){
for(i in 1:length(unique(as_date(clean_UTM_df$datetime)))){

date_i <- unique(as_date(clean_UTM_df$datetime))[i]
date_i_df <- clean_UTM_df[as_date(clean_UTM_df$datetime) == date_i,]
mean_lon <- mean(date_i_df$lon)
mean_lat <- mean(date_i_df$lat)

map = get_map(location = c(lon = mean_lon, lat = mean_lat), zoom=16, maptype="satellite")

g <- ggmap(map)+ 
  geom_path(data=date_i_df, aes(x = lon, y=lat, group= ID, colour= as.factor(name)))+
  geom_point(data=date_i_df, aes(x = lon, y=lat, group=ID, colour=as.factor(name)), size=4)+
  scale_color_manual(values = point_colors) +
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.text=element_text(size=14), 
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))+
  geom_sf(data = dipx, inherit.aes = FALSE, fill = "hotpink", color = "purple", alpha = 0.5) +
  geom_text(data = dipx_centroids, aes(x = X, y = Y, label = newID), size = 3.5, color = "white")

g <- g + transition_reveal(date_i_df$datetime, keep_last = FALSE)+
  ease_aes('linear')+
  labs(title = "Time: {format(as.POSIXct(frame_along, tz = 'UTC'), '%Y-%m-%d %H:%M:%S')}")+
  shadow_wake(wake_length = 0.05, alpha = FALSE)

vid <- gganimate::animate(g, renderer = av_renderer(), height = 1000, width = 1000, fps = 10, duration = 40)

anim_save(filename = paste0(plot_dir, "clean_animations_crowns/V1_movement_", date_i, ".mp4"), animation = vid)
 }
}


#compare distributions of travel speed with the different sampling rates (5 min for adults, 10 min for juveniles)

xs_adult <- clean$xs[coati_ids$age == "Adult",]
ys_adult <- clean$ys[coati_ids$age == "Adult",]

speeds <- matrix(nrow = nrow(xs_adult), ncol = ncol(xs_adult))

for(i in 1:nrow(xs_adult)){

  xs_i <- xs_adult[i,]
  ys_i <- ys_adult[i,]
  
  speed_i <- get_speed(xs_i, ys_i, seconds_per_timestep = 300, t_window = 1)
  speeds[i,] <- speed_i
}

hist(speeds, breaks = 200)
max(speeds, na.rm=TRUE)
median(speeds, na.rm = TRUE)
mean(speeds, na.rm = TRUE)



#now want to look at time spent within dipteryx trees and travel speeds of each group member when at trees - as if the true is not fruiting, they may just travel through, whereas if there's fruit, they are likely going to stay a while to forage
#first interpolating GPS data to have a higher resolution to get more accurate entry and exit times


df <- clean_UTM_df
#df <- df[,-c(3,4)] #remove UTM column
#remove first day as all are every 5 mins but after that only adults are every 5 mins:
df <- df[as.Date(df$datetime) > "2024-12-16",]
colnames(df) <- c("timestamp","tag.local.identifier", "UTM_X", "UTM_Y", "location.long","location.lat","individual.local.identifier")

DATA <- as.telemetry(df, projection = "+proj=utm +zone=17 +north +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

run_code <- F
if(run_code == T){

simulations <- list()

for(ind in names(DATA)) {
  cat("Processing:", ind, "\n")
  
  # Get individual data
  ind_dat <- DATA[[ind]]
  
  # Fit movement model for the individual
  GUESS <- ctmm.guess(ind_dat, interactive = FALSE)
  FIT <- ctmm.select(ind_dat, GUESS, trace = 0)
  
  # Create custom time grid
  t_min <- min(ind_dat$t)
  t_max <- max(ind_dat$t)
  
  time_seq <- seq(from = as.POSIXct(t_min, origin = "1970-01-01", tz = "UTC"),
                  to   = as.POSIXct(t_max, origin = "1970-01-01", tz = "UTC"),
                  by   = "5 min")  # simulate at 2-minute intervals
  
  # Keep only times from 11:00 to 23:00 UTC
  time_seq_filtered <- time_seq[hour(time_seq) >= 11 & hour(time_seq) < 23]
  t_grid_filtered <- as.numeric(time_seq_filtered)
  
  # Simulate using individual model
  SIM <- simulate(ind_dat, FIT, t = t_grid_filtered)
  
  # Add readable time columns
  SIM$time <- as.POSIXct(SIM$t, origin = "1970-01-01", tz = "UTC")
  SIM$date <- as.Date(SIM$time)
  SIM$ID <- ind
  
  simulations[[ind]] <- SIM
}

sim_df <- do.call(rbind, simulations)
save(sim_df, file = paste0(data_dir, "simulated_gps_5min.RData"))

}else if (run_code == F){
  load(paste0(data_dir, "simulated_gps_5min.RData"))
}


dipx_utm <- st_transform(dipx, crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
# Recalculate centroids in UTM
dipx_centroids <- st_centroid(dipx_utm)
dipx_centroids <- cbind(dipx_centroids, st_coordinates(dipx_centroids))  # adds X/Y columns

run_code <- F
if(run_code == T){

plots <- list()

for (ind in unique(sim_df$ID)) {
  sim_i <- sim_df %>% filter(ID == ind)
  obs_i <- df %>% filter(individual.local.identifier == ind)
  
  for (day in unique(sim_i$date)) {
    sim_day <- sim_i %>% filter(date == day)
    obs_day <- obs_i %>% filter(as.Date(timestamp) == day)
    
    day <- as.Date(day)
    
    # Get bounding box of the daily movement (observed + simulated)
    all_xy <- rbind(
      data.frame(x = sim_day$x, y = sim_day$y),
      data.frame(x = obs_day$UTM_X, y = obs_day$UTM_Y)
    )
    
    bbox <- st_as_sfc(st_bbox(c(
      xmin = min(all_xy$x),
      xmax = max(all_xy$x),
      ymin = min(all_xy$y),
      ymax = max(all_xy$y)
    ), crs = st_crs(dipx_utm)))
    
    # Crop Dipteryx trees to just those within this extent
    dipx_crop <- st_intersection(dipx_utm, bbox)
    dipx_centroids_crop <- st_centroid(dipx_crop)
    dipx_centroids_crop <- cbind(dipx_centroids_crop, st_coordinates(dipx_centroids_crop))
    
    p <-  ggplot() +
      geom_path(data = sim_day, aes(x = x, y = y, color = "Simulated"), size = 2, alpha = 0.6) +
      geom_point(data = sim_day, aes(x = x, y = y, color = "Simulated"), size = 0.5, alpha = 1) +
      geom_path(data = obs_day, aes(x = UTM_X, y = UTM_Y, color = "Observed"), size = 1, alpha = 0.8) +
      geom_point(data = obs_day, aes(x = UTM_X, y = UTM_Y, color = "Observed"), size = 0.5, alpha = 0.8) +
      geom_sf(data = dipx_crop, inherit.aes = FALSE, fill = "hotpink", color = "purple", alpha = 0.5) +
      geom_text(data = dipx_centroids_crop, aes(x = X, y = Y, label = newID), size = 1.5, color = "navy") +
      scale_color_manual(name = "", values = c("Observed" = "orangered", "Simulated" = "turquoise2")) +
      labs(title = paste("Observed vs Simulated Movement -", ind, "on", day),
           x = "UTM X", y = "UTM Y") +
      theme_classic() +
      coord_sf() +
      theme(legend.position = "right",
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    
    plots[[paste(ind, day, sep = "_")]] <- p
    
    ggsave(filename = paste0(plot_dir, "/5min_simulation/movement_plot_", ind, "_", day, ".png"),
           plot = p, height = 6, width = 8)
  }
 }
}


























