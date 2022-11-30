#this script is to show where on the map fissions occur

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData' 

library(fields)
library(viridis)
library(ggplot2)
library(ggmap)


#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#load the splits_df made from the split_analysis code in dropbox coatithon 
load("C:/Users/egrout/Dropbox/coatithon/coatithon_code/splits_on_map/splits_df.Rdata")

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
  
#convert to latlon
fission_latlon <- as.data.frame(utm.to.latlon(fission_xy_1, utm.zone = '17',southern_hemisphere=FALSE))
fissions_utm_latlon <- cbind(fission_latlon, fission_xy)
fissions_utm_latlon$event <- as.factor(fissions_utm_latlon$event)

#make map showing where the fission event occured
#NEED TO REGISTER GOOGLE KEY FROM OTHER SCRIPT...
register_google(key="xxx")
has_google_key()
map = get_map(location = c(lon =-79.700003, lat=9.122051), zoom=17, maptype="satellite")
png(file = "C:/Users/egrout/Dropbox/coatithon/results/map_fission_events/fissions_map.png", width = 1000, height = 1000, units = "px")
ggmap(map)+geom_point(data=fissions_utm_latlon, aes(x=V1,y=V2, color=event), size = 4) #+ facet_wrap(~time, ncol = 6)
dev.off()



#-------------------------------------------------------------------------
#get position of individuals before and after split in low res data

#get the index of the moment before there was a split
fission_indx_before <- splits_df$t-1

#find the x and y UTM coords of those splits
fission_xs_before <- xs[,fission_indx_before]
fission_ys_before <- ys[,fission_indx_before]

#make this into a dataframe to reshape it so its not a matrix and that we can join the x and y coords together for the mapping
fission_xs_before <- as.data.frame(fission_xs_before) 
fission_xs_before <- reshape(fission_xs_before, varying=1:ncol(fission_xs_before), v.names="xs", direction ="long", idvar = "ID")

fission_ys_before <- as.data.frame(fission_ys_before)
fission_ys_before <- reshape(fission_ys_before, varying=1:ncol(fission_ys_before), v.names="xs", direction ="long", idvar = "ID")

fission_xy_before <- cbind(fission_xs_before, fission_ys_before)
fission_xy_before <- fission_xy_before[,-c(3,4)]
colnames(fission_xy_before) <- c("event", "xs","ys", "ID")

#quick plot of fission locations 
plot(fission_xy_before$xs, fission_xy_before$ys, col = fission_xy_before$ID) 

#remove id column
fission_xy_1_bef <- fission_xy_before[,-c(1,4)]

#convert to latlon
fission_latlon_before <- as.data.frame(utm.to.latlon(fission_xy_1_bef, utm.zone = '17',southern_hemisphere=FALSE))
fissions_utm_latlon_before <- cbind(fission_latlon_before, fission_xy_before)
fissions_utm_latlon_before$event <- as.factor(fissions_utm_latlon_before$event)

map = get_map(location = c(lon =-79.700003, lat=9.122051), zoom=17, maptype="satellite")
png(file = "C:/Users/egrout/Dropbox/coatithon/results/map_fission_events/before_fissions_map.png", width = 1000, height = 1000, units = "px")
ggmap(map)+geom_point(data=fissions_utm_latlon_before, aes(x=V1,y=V2, color=event), size = 4) 
dev.off()

