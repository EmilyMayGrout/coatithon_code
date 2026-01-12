#Kate's code for looking at coati hotspots during periods of high gps 
library(terra)
library(dplyr)
library(lubridate)
library(sf)

df <- read.csv("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/Nutritional Landscapes_2025_04_21.csv") #movement data from Kate's server, but originally from movebank
#Read in BCI outline 
bciout <- vect("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/bci_outline_new.shp")
dipt <- crowns <- vect("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")

#Get timestamp in panama local time 
df$timestamp = as.POSIXct(strptime(df$timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) #- 3600*5 #Panama time = UTC - 5 h
df$timestamp <- with_tz(df$timestamp, "EST")


#Remove NA fixes
df = subset(df, is.na(df$location.long) != TRUE | is.na(df$location.lat) != TRUE)

#Remove the sloth
df <- df[df$individual.local.identifier != "Jesus",]

#Choose coati or Spider 
taxName <- "Nasua narica"

df <- df[df$individual.taxon.canonical.name == taxName,]


#Choose subset of Days: here I do all of them, you can chose your own 
startday <- force_tz(min(date(df$timestamp)),tz = "EST")
endday <- force_tz(max(date(df$timestamp)),tz = "EST")

#This is insanity, but timezone is a beast  
df <- subset(df, as.Date(df$timestamp, tz = "EST") <= as.Date(endday, tz = "EST") & as.Date(df$timestamp, tz = "EST") >= as.Date(startday, tz = "EST") )

#####No processing done here, but you can and should :) #####
df.sub <- df

#Make it a vect in terra 
dfvect <- vect(df.sub, geom = c("location.long", "location.lat"), crs = "+proj=longlat +datum=WGS84")
dfvect2 <- project(dfvect, crs(bciout))

dfvect2 <- crop(dfvect2, ext(bciout))

dfvect2$Date  <- as.POSIXct(strptime(x=dfvect2$timestamp, format=c("%Y-%m-%d"), tz = 'EST'))
dfvect2$Time <- format(as.POSIXct(dfvect2$timestamp), format = "%H:%M:%S")

#Check to make sure you have right stuff 
table(dfvect2$Date, dfvect2$individual.local.identifier)
#vire 
g <- table(dfvect2$Date, dfvect2$individual.local.identifier)


###If you want just high res periods #####
#eat <- subset(dfvect2, hour( dfvect2$timestamp ) >= 7 & hour(dfvect2$timestamp ) < 11 )                 

tree_crowns_sf <- sf::st_as_sf(crowns)

# Drop empty geometries
tree_crowns_sf <- tree_crowns_sf[!st_is_empty(tree_crowns_sf), ]

# Keep only polygons
tree_crowns_sf <- tree_crowns_sf[
  st_geometry_type(tree_crowns_sf) %in% c("POLYGON", "MULTIPOLYGON"), ]

# Make valid
tree_crowns_sf <- st_make_valid(tree_crowns_sf)

# Drop still-invalid geometries
tree_crowns_sf <- tree_crowns_sf[st_is_valid(tree_crowns_sf), ]

# Drop geometries with less than 4 coordinates (not true polygons)
# Function to count coordinates in each geometry
get_coords_count <- function(geom) {
  coords <- st_coordinates(geom)
  nrow(coords)
}

coord_counts <- sapply(st_geometry(tree_crowns_sf), get_coords_count)
tree_crowns_sf <- tree_crowns_sf[coord_counts >= 4, ]

#make this a spatvect again (sorry!)
crown3 <- vect(tree_crowns_sf)

dfvect2$ID_pt <- seq_len(nrow(dfvect2))

# Extract 
extract_result <- extract(crown3["newID"], dfvect2)

#match as something wonky going on 
point_data <- values(dfvect2)
point_data$ID2 <- seq_len(nrow(dfvect2))

##merge back with pt data and the extracted dipt tree info, usually can cbind 
df_combined <- merge(point_data, extract_result, by.x = "ID2", by.y = "id.y", all.x = TRUE)

head(df_combined)

#remove pts where they didn't intersect with a Dipteryx tree 
dfTreeCrowns <- subset(df_combined, !is.na(df_combined$newID))
table(dfTreeCrowns$newID)

#write your CSV 
write.csv(dfTreeCrowns, paste0("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/", startday, "_", endday, ".csv"))
