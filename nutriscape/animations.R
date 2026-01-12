#making an animation of coati movements in respect to the dipteryx trees

# libraries
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(gganimate)
library(ggmap)
library(ggspatial)
library(ggthemes)
library(sp)
library(sf)
library(fields)
library(viridis)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(plotly)
library(av)
library(cocomo) #from Ari's github
#if need to install cocomo library, run:
#install.packages('devtools’) 
#library(devtools)
#Then run: 
#devtools::install_github('livingingroups/cocomo’) 
#library(cocomo) 


#--------directories-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level0/'
gps_file <- "presidente2025_xy_10min_level0.RData"
id_file <- 'presidente2025_coati_ids.RData'

dipx <- st_read("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#add google map under the points
register_google(key="xxx") #need to get an API for the satellite image, but if you don't want to bother with that, you can make the animations on a different map type
#check you have the google key if using an API
has_google_key()

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds-1)


#-----------------------------------------------------
#visualise movement trajectories 

#function to convert xs and ys matrices to a dataframe for visualising with ggmap
matrix_to_df <- function(xs = NULL, ys = NULL,  ts = NULL, UTM_zone = CRS("+proj=utm +zone=17 +datum=WGS84")){
  
  xs_df <- as.data.frame(xs)
  colnames(xs_df) <- ts
  xs_df$ID <- c(1:16)
  xs_long <- xs_df %>% pivot_longer(!ID, names_to = "datetime", values_to = "UTM_X")
  ys_df <- as.data.frame(ys)
  colnames(ys_df) <- ts
  ys_df$ID <- c(1:16)
  ys_long <- ys_df %>% pivot_longer(!ID, names_to = "datetime", values_to = "UTM_Y")
  UTM_df <- merge(xs_long, ys_long, by = c("datetime", "ID"))
  UTM_df$datetime <- as.POSIXct(UTM_df$datetime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
  
  #add lat and lon to dataframe
  #to do this, need to remove NAs
  UTM_df <- UTM_df[!is.na(UTM_df$UTM_X),]
  
  sp <- SpatialPoints(coords = UTM_df[, c("UTM_X", "UTM_Y")], proj4string = UTM_zone)
  latlon <- spTransform(sp, CRS("+proj=longlat +datum=WGS84"))
  UTM_df <- cbind(UTM_df, unlist(latlon@coords))
  colnames(UTM_df)[c(5,6)] <- c("lon", "lat")
  
  return(UTM_df)
}

pres_df <- matrix_to_df(xs, ys, ts)


#plot all raw data to see outliers
plot(pres_df$UTM_X, pres_df$UTM_Y, type = "p", asp = 1, cex = 0.1)


#clean the data to remove GPS outliers
breaks <- (which(diff(as_date(ts)) == 1)) + 1 #adding 1 as it gets the last each date and we need the index of the first time point for each day

clean <- preprocess_gps_level0_to_level1(xs = xs, ys = ys, timestamps = ts, ids = coati_ids, breaks = breaks, max_isolated_point_dist = 85, max_dist_percentile = 0.98, verbose = F)


#use function made above to convert matrix to dataframe for plotting the clean data
clean_UTM_df <- matrix_to_df(xs = clean$xs, ys = clean$ys, ts = clean$timestamps)

#adding ID column to combine the coati ids with the clean dataframe
coati_ids$ID <- c(1:16)

clean_UTM_df <- clean_UTM_df %>%
  left_join(coati_ids[, c("ID", "name")], by = "ID")

point_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928', "black", "white", "grey", "yellow", "darkblue")


map = get_map(location = c(lon = mean(clean_UTM_df$lon), lat = mean(clean_UTM_df$lat)), zoom=16, maptype="satellite")

#plot all the data on satellite image to see the hotspots 
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

#transform the dipteryx shapefile to the same CRS
dipx <- st_transform(dipx, crs = 4326)

#get bounding box of shapefile
bbox <- st_bbox(dipx)

#make a base map
basemap <- get_map(
  location = c(lon = mean(clean_UTM_df$lon), 
               lat = mean(clean_UTM_df$lat)),
  zoom = 16,   # adjust zoom level as needed
  maptype = "satellite"  # or "satellite", "roadmap", etc.
)


#fix invalid geometries
dipx <- st_make_valid(dipx)

#calculate centroids for where the tree IDs will be written
dipx_centroids <- st_centroid(dipx)

dipx_centroids <- dipx_centroids %>%
  st_transform(4326) %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  )


#remove Ardern from visualisations as she's on a lower res
clean_UTM_df <- clean_UTM_df[clean_UTM_df$name != "Ardern",]

#this is running a for loop through each day so takes hours to run, if you just want to make one (which should be around 15 mins, just run the code without the for loop and choose an i value)
run_code <- F

if(run_code == T){
  for(i in 1:length(unique(as_date(clean_UTM_df$datetime)))){
    
    #i <- 1
    
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


