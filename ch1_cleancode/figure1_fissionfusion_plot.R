#this code is to make the before and after fission plot for figure 1 (panel b and c)
#here is also the code for the full group trajectories for 1 day

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/'
gps_file <- "galaxy_xy_10min_level1.RData"
id_file <- 'galaxy_coati_ids.RData' 

#load in libraries
library(fields)
library(viridis)
library(sf)
library(purrr)
library(ggplot2)
library(sf)
library(ggmap)
library(dplyr)
library(lubridate)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#load in splits df from fission_fusion_galaxy_V1 to get the split times, so before and after the event can be extracted for the ff plot
load(file = "C:/Users/egrout/Dropbox/coatithon_notgithub/Galaxy_fission_fusion/splits_df.Rdata")

#-----MAIN------
n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, 50)

#get the index of when there was a split
fission_indx <- splits_df$t

#find the x and y UTM coords of those splits
fission_xs <- xs[,fission_indx]
fission_ys <- ys[,fission_indx]

#make this into a dataframe to reshape it so its not a matrix and that we can join the x and y coords together for the mapping
fission_xs <- as.data.frame(fission_xs) 
fission_xs <- reshape(fission_xs, varying=1:ncol(fission_xs), v.names="xs", direction ="long", idvar = "ID")

fission_ys <- as.data.frame(fission_ys)
fission_ys <- reshape(fission_ys, varying=1:ncol(fission_ys), v.names="xs", direction ="long", idvar = "ID")

fission_xy <- cbind(fission_xs, fission_ys)
fission_xy <- fission_xy[,-c(3,4)]
colnames(fission_xy) <- c("event", "xs","ys", "ID")

#quick plot of fission locations 
plot(fission_xy$xs, fission_xy$ys, col = fission_xy$ID) 

#remove id column
fission_xy_1 <- fission_xy[,-c(1,4)]

# #convert to latlon

#convert NA's to zero for sf functions to work
fission_xy_1$xs[is.na(fission_xy_1$xs)] <- 0
fission_xy_1$ys[is.na(fission_xy_1$ys)] <- 0

fission_utm <- st_as_sf(x=fission_xy_1, coords=c("xs", "ys"), crs="+proj=utm +zone=17 +north +datum=WGS84 +units=m")
#convert to UTM
fission_latlon <- st_transform(fission_utm, crs= "+proj=longlat +datum=WGS84") #convert UTM to lat/long - coordinates are stored in a geometry column of class 'sfc_POINT'

#store lat and long in the dataframe with the id and utm coords
lon <- unlist(purrr::map(fission_latlon$geometry,1)) #longitude is X
lat <- unlist(purrr::map(fission_latlon$geometry,2)) #latitude is Y

fission_utm_latlon <- cbind(fission_xy, lon, lat)

#change the rows to 0 if lat or lon values are 0
fission_utm_latlon$lon[fission_utm_latlon$lat == 0] <- 0
fission_utm_latlon[fission_utm_latlon == 0] <-  NA

fission_utm_latlon$event <- as.factor(fission_utm_latlon$event)


#-------------------------------------------------------
#get the xs and ys for the times of the splits
split_xs <- xs[, splits_df$t]
before_split_xs <- xs[,(splits_df$t)-1]
after_split_xs <- xs[,(splits_df$t)+1]
split_ys <- ys[, splits_df$t]
before_split_ys <- ys[,(splits_df$t)-1]
after_split_ys <- ys[,(splits_df$t)+1]

#-----------------------------------------------------------------------------------
#combine the xs and ys of before into the merged_df and convert them to latlon
xs_bef <- as.data.frame(before_split_xs)
ys_bef <- as.data.frame(before_split_ys)
xs_aft <- as.data.frame(after_split_xs)
ys_aft <- as.data.frame(after_split_ys)

# Get the row names (individuals) and add them as a new column
xs_bef$id <- rownames(xs_bef)
ys_bef$id <- rownames(ys_bef)
xs_aft$id <- rownames(xs_aft)
ys_aft$id <- rownames(ys_aft)

# Reshape the data from wide to long format
xs_bef_long <- reshape::melt(xs_bef)
ys_bef_long <- reshape::melt(ys_bef)
xs_aft_long <- reshape::melt(xs_aft)
ys_aft_long <- reshape::melt(ys_aft)

colnames(xs_bef_long) <- c("id", "event", "xs_bef")
colnames(ys_bef_long) <- c("id", "event", "ys_bef")
colnames(xs_aft_long) <- c("id", "event", "xs_aft")
colnames(ys_aft_long) <- c("id", "event", "ys_aft")

#remove "V" from event column
xs_bef_long$event <- gsub("V", "", xs_bef_long$event)
ys_bef_long$event <- gsub("V", "", ys_bef_long$event)
xs_aft_long$event <- gsub("V", "", xs_aft_long$event)
ys_aft_long$event <- gsub("V", "", ys_aft_long$event)

split_df_all <- cbind(fission_utm_latlon, xs_bef_long, ys_bef_long, xs_aft_long, ys_aft_long)

split_df_all <- split_df_all[, -c(7,8,10,11,13,14,16,17)]

#get latlon for splits_bef and split_aft
#convert NA's to zero for sf functions to work

split_df_all$xs_bef[is.na(split_df_all$xs_bef)] <- 0
split_df_all$ys_bef[is.na(split_df_all$ys_bef)] <- 0

split_utm <- st_as_sf(x=split_df_all, coords=c("xs_bef", "ys_bef"), crs="+proj=utm +zone=17 +north +datum=WGS84 +units=m")
#convert to UTM
split_latlon <- st_transform(split_utm, crs= "+proj=longlat +datum=WGS84") #convert UTM to lat/long - coordinates are stored in a geometry column of class 'sfc_POINT'

#store lat and long in the dataframe with the id and utm coords
lon_bef <- unlist(map(split_latlon$geometry,1)) #longitude is X
lat_bef <- unlist(map(split_latlon$geometry,2)) #latitude is Y

split_df_all <- cbind(split_df_all, lon_bef, lat_bef)

#change the rows to 0 if lat or lon values are 0
split_df_all$lon_bef[split_df_all$lat_bef == 0] <- 0
split_df_all[split_df_all == 0] <-  NA

#also do for aft xs and ys
split_df_all$xs_aft[is.na(split_df_all$xs_aft)] <- 0
split_df_all$ys_aft[is.na(split_df_all$ys_aft)] <- 0

split_utm <- st_as_sf(x=split_df_all, coords=c("xs_aft", "ys_aft"), crs="+proj=utm +zone=17 +north +datum=WGS84 +units=m")
#convert to UTM
split_latlon <- st_transform(split_utm, crs= "+proj=longlat +datum=WGS84") #convert UTM to lat/long - coordinates are stored in a geometry column of class 'sfc_POINT'

#store lat and long in the dataframe with the id and utm coords
lon_aft <- unlist(map(split_latlon$geometry,1)) #longitude is X
lat_aft <- unlist(map(split_latlon$geometry,2)) #latitude is Y

split_df_all <- cbind(split_df_all, lon_aft, lat_aft)

#change the rows to 0 if lat or lon values are 0
split_df_all$lon_aft[split_df_all$lat_aft == 0] <- 0
split_df_all[split_df_all == 0] <-  NA

#split with dbscan
event_i <- split_df_all[split_df_all$event == 25,]

#getting max and min lat lons for a buffer for the map
min_lat <- (min(c(event_i$lat, event_i$lat_bef,  event_i$lat_aft), na.rm = TRUE) - 0.0005)
max_lat <- (max(c(event_i$lat, event_i$lat_bef,  event_i$lat_aft), na.rm = TRUE) + 0.0005)
min_lon <- (min(c(event_i$lon, event_i$lon_bef, event_i$lon_aft), na.rm = TRUE) - 0.0005)
max_lon <- (max(c(event_i$lon, event_i$lon_bef, event_i$lon_aft), na.rm = TRUE) + 0.0005)

#for making a scale bar
sites.data = data.frame(lon = c(min_lon, max_lon),
                        lat = c(min_lat, max_lat))

#make your map
map.base <- get_map(location = c(lon = mean(sites.data$lon),
                                 lat = mean(sites.data$lat)), zoom = 18)

#ggmap(map.base)
bb <- attr(map.base,"bb")
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                  c(sbar$lon.end,sbar$lat.end))

scalebar.length <- 0.05
sbar$lon.end <- sbar$lon.start + ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000

#------------------------------------------------------------------------------------------
#FIG1: make plot for group split with DBSCAN for each individual

lon_time <- event_i$lon
lat_time <- event_i$lat

#for event 25 (for paper, used size 130 for aft split and 120 for before (due to change in plot margins))

gg <- #ggmap(map)+
  #if don't want on map background, change to ggplot and add panel info in theme, changed size of circles to 115 from 75 
  ggplot(event_i, aes(x=lon_time,y=lat_time))+
  geom_point(data=event_i, aes(x=lon_time,y=lat_time), color = "white", size = 120, alpha = 0.5, shape = 1)+ #outline circle
  geom_point(data=event_i, aes(x=lon_time,y=lat_time), color = "white", size = 120, alpha = 0.1)+ #outline circle
  geom_point(data=event_i, aes(x=lon_time,y=lat_time), color = "white", size = 9, alpha = 0.8)+
  scale_x_continuous(limits = c(min_lon, max_lon))+
  scale_y_continuous(limits = c(min_lat, max_lat))+ 
  xlab(" ") +  # Add X-axis label
  ylab(" ") +  # Add Y-axis label
  ggtitle("(b)")+
  annotate("text", x=max_lon - 0.0001, y=min_lat, label= "t", color = "yellow", size = 20)+ #for aft, do max_lon - 0.0001, for dur, only max_lon
  #adding a scale bar
  geom_segment(data = sbar,
               aes(x = lon.start + 0.0018,
                   xend = lon.end + 0.0018,
                   y = lat.start+0.0022,
                   yend = lat.end+0.0022),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open"), color = "white", lwd=3) +
  geom_text(data = sbar,
            aes(x = (lon.start +0.0018 + lon.end +0.0018)/2,
                y = lat.start+ 0.0023 + 0.001*(bb$ur.lat - bb$ll.lat),
                label = '50 m'),
            hjust = 0.5,
            vjust = 0.7,
            size = 12, color = "white") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "turquoise4",colour = "turquoise4", linetype = "solid"),
        plot.background = element_rect(fill = "turquoise4"),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=60, color = "white", hjust = 0.03, vjust = -1))

gg

ggsave(filename = paste0(plot_dir, 'ggmap_split_bef', 25, '.png'), plot = gg, width = 11, height = 10, dpi = 300)

#------------------------------------------------------------------------------------------
#make plot of groups trajectories using the lat/lon data from preprocess_gps_lowres_galaxy.R script
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_latlon_10min_level0.RData")

#put into a dataframe format
df_i <- data.frame(id = 1:ncol(lats), lat = 1:ncol(lats), lon = 1:ncol(lats), ts = 1:ncol(lats))

df <- data.frame(id = NA, lat = NA, lon =NA, ts = NA)

for (i in 1:nrow(lats)){
  
  #make the row i of the lats into a dataframe to add to the lat_df dataframe
  lat_i <- lats[i,]
  lon_i <- lons[i,]
  #add it to a dataframe with the id and time
  df_i$lat <- lat_i
  df_i$lon <- lon_i
  df_i$ts <- ts
  df_i$id <- i
  
  #rbind ind i to the lat_df dataframe (which will be used for the map plots)
  df <- rbind(df, df_i)
  
  #clean the df_i dataframe for next i 
  df_i <- data.frame(id = 1:ncol(lats), lat = 1:ncol(lats), lon = 1:ncol(lats), ts = 1:ncol(lats))
  
}

#remove rows with NA
df <- na.omit(df)

#filter data frame to one date
df$ts <- as.POSIXct(df$ts)
df_subset <- subset(df, ts > as.POSIXct("2021-12-28 14:00:00") & ts < as.POSIXct("2021-12-28 23:00:00"))
df_subset <- subset(df_subset, id != 5) 

#plot the lat and lons for all points
map = get_map(location = c(lon =-79.700063, lat=9.123001), zoom=15, maptype = "satellite")
#getting max and min lat lons for a buffer for the map
min_lat <- (min(df_subset$lat, na.rm = TRUE) - 0.0001)
max_lat <- (max(df_subset$lat, na.rm = TRUE) + 0.0000) #changed from 0.0020
min_lon <- (min(df_subset$lon, na.rm = TRUE) - 0.0005)
max_lon <- (max(df_subset$lon, na.rm = TRUE) + 0.0005)
#for making a scale bar
sites.data = data.frame(lon = c(min_lon, max_lon),
                        lat = c(min_lat, max_lat))

map.base <- get_map(location = c(lon = mean(sites.data$lon),
                                 lat = mean(sites.data$lat)), zoom = 18)

#ggmap(map.base)
bb <- attr(map.base,"bb")
sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                   lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                   lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                   lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))

sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                  c(sbar$lon.end,sbar$lat.end))

scalebar.length <- 0.05
sbar$lon.end <- sbar$lon.start + ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000

df_subset$id <- as.factor(df_subset$id)
# Manually specify colors for id levels
#id_colours <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#990000", "#009900", "#000099", "#999900", "#990099")


##for this plot to work, you need to register a google key
#register_google(key="xxx")
has_google_key()

df_subset$group <- 1

row.names(df_subset) <- 1:nrow(df_subset)

#plot for group trajectories
gg <- ggmap(map) +
  geom_path(data = unique(df_subset), aes(x = lon, y = lat, group = id), color = "white", size = 0.2, alpha = 1) +
  geom_point(data = df_subset, aes(x = lon, y = lat, color = ts), size = 3.5, alpha = 0.6) +
  scale_colour_gradient(low = "purple", high = "yellow")+
  #data = df_subset, aes(x = lon, y = lat, color = as.factor(id)), size = 2, alpha = 0.6) #scale_color_manual(values = id_colours)+
  
  scale_x_continuous(limits = c(min_lon, max_lon))+
  scale_y_continuous(limits = c(min_lat, max_lat))+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  #adding a scale bar
  geom_segment(data = sbar,
               aes(x = lon.start + 0.0038,
                   xend = lon.end + 0.0038,
                   y = lat.start+0.0015,
                   yend = lat.end+0.0015),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open"), color = "white") +
  geom_text(data = sbar,
            aes(x = (lon.start +0.0038 + lon.end +0.0038)/2,
                y = lat.start+ 0.0017 + 0.001*(bb$ur.lat - bb$ll.lat),
                label = '50m'),
            hjust = 0.5,
            vjust = 0.7,
            size = 8, color = "white") +
  NULL
gg


ggsave(filename = paste0(plot_dir, 'ggmap_traj_28', '.png'), plot = gg, width = 14, height = 10, dpi = 300)

#--------------------------------------------------------------------------

#add column with subgroup id
df_subset <- df_subset %>%
  mutate(subgroup_id = ifelse(id %in% c(3,6,7,8,9), 1, ifelse(id %in% c(1,2,4,5,10,11), 2, NA)))

#trying different plot with colours representing the subgroups
gg <- ggmap(map) +
  geom_path(data = df_subset, aes(x = lon, y = lat, group = id), color = "white", size = 1, alpha = 0.5)+
  geom_point(data = df_subset, aes(x = lon, y = lat, color = subgroup_id), size = 3, alpha = 0.6)+
  scale_color_gradient(low="mediumpurple1", high="cadetblue1") +
 # geom_point(data = df_subset, aes(x = lon, y = lat), color = "white", shape = 1, size = 2, alpha = 1)+
  scale_x_continuous(limits = c(min_lon, max_lon))+
  scale_y_continuous(limits = c(min_lat, max_lat))+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  #adding a scale bar
  geom_segment(data = sbar,
               aes(x = lon.start + 0.0038,
                   xend = lon.end + 0.0038,
                   y = lat.start+0.0015,
                   yend = lat.end+0.0015),
               arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                           ends = "both", type = "open"), color = "white") +
  geom_text(data = sbar,
            aes(x = (lon.start +0.0038 + lon.end +0.0038)/2,
                y = lat.start+ 0.0017 + 0.001*(bb$ur.lat - bb$ll.lat),
                label = '50 m'),
            hjust = 0.5,
            vjust = 0.7,
            size = 8, color = "white") +
  #geom_text(data = sbar,
  #          aes(x = -79.7052,
  #              y = 9.1235,
  #              label = '(a)'), color = "white", size = 8)+
  
  NULL

gg

ggsave(filename = paste0(plot_dir, 'ggmap_traj_28_col_cut2', '.png'), plot = gg, width = 13, height = 6, dpi = 300)






