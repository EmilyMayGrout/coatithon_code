#this script is reading in the ACC and taking a first look at the data
#goal is to extract vedba and look at how this changes when coatis are at hotspots - is there variation between group members when entering a food patch to leaving?

library(dplyr)
library(ggplot2)
library(lubridate)
library(mclust)
library(dbscan)
library(sf)
library(patchwork)
library(cowplot)
library(hms)

plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level1/'

#load in all movebank data - 14.26 end 14.30
#all <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_all.csv")
#filter to coati data
#taxName <- "Nasua narica"
#coati_all <- all[all$individual.taxon.canonical.name == taxName,]
#save(coati_all, file = "C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_all.RData")

# load("C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_all.RData")
#filter to just ACC data
#coati_acc <- coati_all[coati_all$sensor.type == "acceleration",]

# coati_acc_vedba <- coati_acc %>%
#   rowwise() %>%
#   mutate(vedba = list({
#     acc_values <- as.numeric(unlist(strsplit(eobs.accelerations.raw, " ")))
#     acc_matrix <- matrix(acc_values, ncol = 3, byrow = TRUE)
#     acc_static <- apply(acc_matrix, 2, function(col) zoo::rollmean(col, k = 20, fill = NA, align = "center"))
#     acc_dynamic <- acc_matrix - acc_static
#     vedba_series <- sqrt(rowSums(acc_dynamic^2, na.rm = TRUE))
#     mean(vedba_series, na.rm = TRUE)
#   }))
#
# save(coati_acc_vedba, file = "C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_vedba.RData")

load("C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_vedba.RData")

#2 min GPS simulated from the 5 and 10 minute sampling in fist_look_at_data script
load("C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/simulated_gps_2min.RData")
#dipteryx shapefile
dipx <- st_read("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")

#filter to just columns we need and rename columns to match to GPS df
columns_to_keep <- c("timestamp", "vedba", "individual.local.identifier", "eobs.accelerations.raw")
coati_vedba <- coati_acc_vedba
coati_vedba <- coati_vedba[ , columns_to_keep]
coati_vedba$vedba <- as.numeric(coati_vedba$vedba)
#round ACC to nearest 10 second to match with GPS (only moves the time 1s back)
coati_vedba$timestamp <- as.POSIXct(coati_vedba$timestamp, origin = "1970-01-01", tz = "UTC")
coati_vedba$timestamp <- floor_date(coati_vedba$timestamp, unit="60 second")
names(coati_vedba)[names(coati_vedba) == 'individual.local.identifier'] <- 'ID'

hist(coati_vedba$vedba, breaks = 100)
hist(log(coati_vedba$vedba), breaks = 100)

ggplot(coati_vedba, aes(x = log(vedba)))+
  geom_histogram(bins = 100)+
  facet_wrap(~ID)+
  geom_vline(xintercept = log(67))+
  #xlim(0,300)+
  #ylim(0, 2500)+
  theme_classic()

#why does Einstein have less data?
#Einstein <- coati_vedba[coati_vedba$individual.local.identifier == "Einstein",]
#his ACC stopped after the 26th of December
#max(Einstein$timestamp)

#look at vedba distributions for day period only:
coati_vedba$hour <- hour(coati_vedba$timestamp)

vedba_daytime <- coati_vedba[which((coati_vedba$hour <= 23) & (coati_vedba$hour >= 11)),]
vedba_nighttime <- coati_vedba[which((coati_vedba$hour <= 10) & (coati_vedba$hour >= 1)),]

ggplot(vedba_daytime, aes(x= log(vedba)))+
  geom_histogram(bins = 100)+
  facet_wrap(~ID)

#looking at just daytime vedba to see how this changes the thresholds
#coati_vedba <- vedba_daytime

#----------------------------------------------------------
#extract the distribution of the Z axis per individual to look at when coatis have entered a fruit tree, did they climb the tree

# acc_list <- strsplit(coati_vedba$eobs.accelerations.raw, " ")
# 
# # Compute mean and sd for Z (3rd column)
# z_stats <- lapply(acc_list, function(x) {
#   vals <- as.numeric(x)
#   z_vals <- matrix(vals, ncol = 3, byrow = TRUE)[, 3]
#   c(mean = mean(z_vals), sd = sd(z_vals))
# })
# 
# # Bind back into a dataframe
# z_stats_df <- do.call(rbind, z_stats)
# 
# coati_vedba$mean_z <- z_stats_df[, "mean"]
# coati_vedba$sd_z   <- z_stats_df[, "sd"]
# 
# coati_vedba_z <- coati_vedba
# 
# save(coati_vedba_z, file = "C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_vedba_zaxis.RData")

load("C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes_vedba_zaxis.RData")

hist(coati_vedba_z$mean_z, breaks= 100)

#--------------------------------------------------------------

#merge gps with the acc dataframe
sim_df$timestamp <- as.POSIXct(sim_df$t, origin = "1970-01-01", tz = "UTC")
#filter to just columns needed
columns_to_keep <- c("timestamp","x", "y","ID")
coati_gps <- sim_df[ , columns_to_keep]
coati_gps$dates <- as.Date(coati_gps$timestamp)

gps_summary <- coati_gps %>%
  group_by(ID) %>%
  summarise(
    start_date = min(dates),
    end_date = max(dates),
    n_days = n_distinct(dates)
  ) %>%
  arrange(desc(n_days))  # optional: sort by number of days

# Plot number of days recorded per individual
ggplot(gps_summary, aes(x = reorder(ID, -n_days), y = n_days)) +
  geom_col(fill = "steelblue") +
  labs(title = "Number of Days with GPS Data per Individual",
       x = "Individual",
       y = "Number of Days") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

columns_to_keep <- c("timestamp","x", "y","ID")
coati_gps <- coati_gps[ , columns_to_keep]
coati_gps$local_timestamp <- coati_gps$timestamp - as.difftime(5, units = "hours")

#so the IDs are the same for merging acc and gps data
coati_vedba$ID[coati_vedba$ID == "Da Vinci"] <- "DaVinci"

#look at proportion of time each individual is active
#first need to find cut-off for active vs resting - used a gausian mixture model for each individual
#fit a mixture of two normal distributions to the VeDBA values, and find the intersection (i.e., the value where the two distributions cross)
vedba_thresholds <- coati_vedba %>%
  group_by(ID) %>%
  mutate(
    threshold = mean(mclust::Mclust(vedba, G = 2)$parameters$mean), 
    activity_level = ifelse(vedba > threshold, "high", "low")
  )

#how to get 2 thresholds for low, medium, and high activity levels as the current thresholds for the daytime data are in the middle of the larger peaks (where there isn't a clear cut-off point)
#might be because general higher activity during the day, so threshold needs more low activity data - could use the thresholds calculated with all day and night data to find the general activity pattern thresholds instead, and then subsample to the day period?

unique(vedba_thresholds$threshold)

#plot the data with each individuals cutoff
ggplot(vedba_thresholds, aes(x = vedba)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  geom_vline(aes(xintercept = threshold), color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ ID, scales = "free_y") +
  xlim(0,500)+
  theme_minimal() +
  labs(
    title = "VeDBA Distributions with Threshold per Individual",
    x = "VeDBA",
    y = "Density")

merged_data <- full_join(vedba_thresholds, coati_gps, by = c("timestamp", "ID"))

#get the proportion of time each individual is in low vs high vedba
prop_vedba <- as.data.frame(table(vedba_thresholds$ID, vedba_thresholds$activity_level))
colnames(prop_vedba) <- c("ID", "activity_level", "count")
total_counts <- aggregate(count ~ ID, data = prop_vedba, FUN = sum)
colnames(total_counts)[2] <- "total"
activity_proportions <- merge(prop_vedba, total_counts, by = "ID")
activity_proportions$proportion <- activity_proportions$count / activity_proportions$total

# Create a matrix of proportions with individuals as rows and activity levels as columns
prop_matrix <- reshape(activity_proportions[, c("ID", "activity_level", "proportion")], timevar = "activity_level", idvar = "ID", direction = "wide")

# Set row names to individual IDs
rownames(prop_matrix) <- prop_matrix$ID
prop_matrix <- prop_matrix[, -1]  # remove ID column

# Rename columns for clarity
colnames(prop_matrix) <- c("High", "Low")

#reorder for plotting low below high 
prop_matrix <- prop_matrix[,c(2,1)]

# Transpose for stacked barplot
prop_matrix_t <- t(as.matrix(prop_matrix))

# Create stacked bar plot
barplot(prop_matrix_t,
        beside = FALSE,
        col = c("skyblue", "tomato"),
        legend = TRUE,
        args.legend = list(title = "Activity Level", x = "topright"),
        ylab = "Proportion of Time",
        xlab = "",
        main = "Proportion of Time in Low vs High Activity",  las = 2)

#how does this change for each hour of the day?
vedba_thresholds$local_timestamp <- vedba_thresholds$timestamp - as.difftime(5, units = "hours")
vedba_thresholds$hour <- hour(vedba_thresholds$local_timestamp) #for local time for plot

activity_by_hour <- vedba_thresholds %>%
  group_by(ID, hour, activity_level) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ID, hour) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()


ggplot(activity_by_hour, aes(x = hour, y = proportion, fill = activity_level)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ ID, ncol = 4) +  # Adjust columns as needed
  scale_fill_manual(values = c("low" = "skyblue3", "high" = "yellow")) +
  labs(
    title = "Hourly Activity Levels by Individual",
    x = "Hour of Day",
    y = "Proportion of Time",
    fill = "Activity Level"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(hjust = 1)
  )
#next step is to isolate times when individuals are in a hotspot and look at their activity patterns. 
#as the fruits of a dipteryx are sometimes not directly below the tree, should use the hotspots the coatis spend lots of time in and then later see whether they were at or near a dipteryx tree or perhaps it was a different fruit tree, or a rest site. Should also see where the sleep trees are and look at the groups sleeping behaviour throughout nights 

point_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928', "black", "white", "grey", "yellow", "darkblue")

#look at the dipteryx tree crowns in relation to the hotspots

# Get bounding box of shapefile
bbox <- st_bbox(dipx)

# First, fix invalid geometries
dipx <- st_make_valid(dipx)

# Then calculate centroids
dipx_centroids <- st_centroid(dipx)

dipx_centroids <- dipx_centroids %>%
  mutate(
    X = st_coordinates(.)[, 1],
    Y = st_coordinates(.)[, 2]
  ) 

run_code <- F

if(run_code == T){
for (i in 1:length(unique(as.Date(coati_gps$local_timestamp)))){

#look at one day to isolate hotspots
date <- unique(as.Date(coati_gps$local_timestamp))[i]
df_cut <- coati_gps[as.Date(coati_gps$local_timestamp) == date,]

#find sleep site locations
first_time <- min(df_cut$local_timestamp)
morning_sleep_site <- df_cut[df_cut$local_timestamp == first_time,]
morning_sleep_site <- morning_sleep_site %>% group_by(ID) %>%
summarise(mean(x), mean(y))
colnames(morning_sleep_site) <- c("ID", "mean_x", "mean_y")

#find clusters where there are at least 80 points within a 20m radius
db <- dbscan(df_cut[, c("x", "y")], eps = 20, minPts = 80) 

df_cut$cluster <- db$cluster
hotspot_means <- df_cut %>%
  filter(cluster != 0) %>%
  group_by(cluster) %>%
  summarise(mean_x = mean(x), mean_y = mean(y))

x_range <- range(df_cut$x, na.rm = TRUE)
y_range <- range(df_cut$y, na.rm = TRUE)

g <- ggplot(df_cut, aes(x = x, y = y, group=ID, colour= as.factor(ID)))+
  geom_path(alpha = 0.3) +
  scale_color_manual(values = point_colors) +
  coord_fixed(ratio = 1) +
  labs(title = paste0("Hotspots (dbscan, R = 20m, minPts = 80) for ", date))+
  theme_classic() + 
  geom_point(data = morning_sleep_site, aes(x=mean_x, y=mean_y), inherit.aes = FALSE, colour = "slateblue4", shape = 16)+
  #theme(legend.position="none")+ 
  geom_sf(data = dipx, inherit.aes = FALSE, fill = "pink", alpha = 0.5, color = "purple4") +
 geom_point(data = hotspot_means, aes(x = mean_x, y = mean_y), color = "cyan", shape = 4, size = 4, inherit.aes = FALSE, stroke = 1.5)+ 
  geom_text(data = dipx_centroids, inherit.aes = FALSE,
            aes(x = X, y = Y, label = newID),
            size = 1.5, color = "red3")+
  xlim(x_range) +
  ylim(y_range)
g

ggsave(filename = paste0(plot_dir, 'hotspots/hotspots20m_80minpts/hotspots_', date, '.png'), plot = g, width = 8, height = 8, dpi = 300)

 }
}

#look at one day of the gps and plotting the vedba 
merged_data <- merge(coati_gps, vedba_thresholds[, c("local_timestamp", "ID", "activity_level")], by = c("local_timestamp", "ID"))


if(run_code == T){
for (i in 1:length(unique(as.Date(merged_data$local_timestamp)))){
  
  date <- unique(as.Date(merged_data$local_timestamp))[i]
  cut_data <- merged_data[as.Date(merged_data$local_timestamp) == date,]
  x_range <- range(cut_data$x, na.rm = TRUE) 
  y_range <- range(cut_data$y, na.rm = TRUE) 
  first_time <- min(cut_data$local_timestamp)
  morning_sleep_site <- cut_data[cut_data$local_timestamp == first_time,]
  morning_sleep_site <- morning_sleep_site %>% group_by(ID) %>%
    summarise(mean(x), mean(y))
  colnames(morning_sleep_site) <- c("ID", "mean_x", "mean_y")
  rest_spots <- cut_data[cut_data$activity_level == "low",]
  
  #get the location hotspots
  db <- dbscan(cut_data[, c("x", "y")], eps = 20, minPts = 100) 
  cut_data$cluster <- db$cluster
  hotspot_means <- cut_data %>%
    filter(cluster != 0) %>%
    group_by(cluster) %>%
    summarise(mean_x = mean(x), mean_y = mean(y))
  
g <- ggplot(cut_data, aes(x = x, y = y, color = activity_level)) +
  geom_path(data = cut_data, aes(x = x, y = y, group = ID), 
            color = "yellow4", alpha = 0.4) +
  geom_point(alpha = 0.2) +
 
  geom_point(data = cut_data[cut_data$activity_level == "high", ], 
             aes(x = x, y = y), color = "tomato", alpha = 0.2) +
  geom_point(data = cut_data[cut_data$activity_level == "low", ], 
             aes(x = x, y = y), color = "skyblue", alpha = 0.6) + 
  
  coord_fixed() +
  labs(title = paste("Activity Locations on", date), 
       color = "Activity Level") +
  geom_sf(data = dipx, inherit.aes = FALSE, fill = "pink", 
          alpha = 0.5, color = "purple3") +
  geom_point(data = hotspot_means, aes(x = mean_x, y = mean_y), color = "cyan", shape = 4, size = 4, inherit.aes = FALSE, stroke = 2)+ 
  geom_point(data = morning_sleep_site, aes(x=mean_x, y=mean_y), inherit.aes = FALSE, 
             colour = "slateblue4", shape = 4)+
  geom_text(data = dipx_centroids, inherit.aes = FALSE,
            aes(x = X, y = Y, label = newID),
            size = 1.5, color = "violetred4") +
  xlim(x_range) +
  ylim(y_range) +
  theme_classic()

g

ggsave(filename = paste0(plot_dir, 'hotspots/vedba/vedba_', date, '.png'), plot = g, width = 9, height = 8, dpi = 300)

}
}



#-------------------------------------------

#plotting the vedba values when coatis are in each cluster per day


# Assign consistent colors
unique_ids <- unique(coati_vedba$ID)
custom_colors <- c("darkorange", "darkgreen", "steelblue", "orchid", "brown", 
                   "dodgerblue", "deeppink", "gold", "slateblue", "firebrick", 
                   "seagreen", "navy", "mediumpurple", "coral", "grey40")
id_colors <- setNames(custom_colors[1:length(unique_ids)], unique_ids)

# Get high activity threshold
thresh <- mean(unique(vedba_thresholds$threshold))

# Loop over each date
for (i in 1:length(unique(as.Date(coati_gps$local_timestamp)))) {
  
  date <- unique(as.Date(coati_gps$local_timestamp))[i]
  cut_data <- coati_gps[as.Date(coati_gps$local_timestamp) == date,]
  x_range <- range(cut_data$x, na.rm = TRUE) 
  y_range <- range(cut_data$y, na.rm = TRUE) 
  
  morning_sleep_site <- cut_data %>%
    filter(local_timestamp == min(local_timestamp)) %>%
    group_by(ID) %>%
    summarise(mean_x = mean(x), mean_y = mean(y))
  
  db <- dbscan(cut_data[, c("x", "y")], eps = 20, minPts = 80) 
  cut_data$cluster <- db$cluster
  
  hotspot_means <- cut_data %>%
    filter(cluster != 0) %>%
    group_by(cluster) %>%
    summarise(mean_x = mean(x), mean_y = mean(y))
  
  hotspot_visits <- cut_data %>%
    filter(cluster != 0) %>%
    group_by(ID, cluster) %>%
    summarise(
      entry_time = min(local_timestamp),
      exit_time = max(local_timestamp),
      duration_min = as.numeric(difftime(exit_time, entry_time, units = "mins")),
      .groups = "drop"
    )
  
  mean_durations <- hotspot_visits %>%
    group_by(cluster) %>%
    summarise(mean_duration = round(mean(duration_min), 1))
  
  hotspot_labels <- hotspot_means %>%
    left_join(mean_durations, by = "cluster") %>%
    mutate(label = paste0("Cluster: ", cluster, "\nMean: ", mean_duration, " min"))
  
  map_plot  <- ggplot(cut_data, aes(x = x, y = y, color = ID)) +
    geom_path(aes(group = ID), color = "grey", alpha = 0.4) +
    geom_point(alpha = 0.2) +
    scale_color_manual(values = id_colors) +
    geom_point(data = hotspot_means, aes(x = mean_x, y = mean_y), 
               color = "cyan", shape = 4, size = 10, inherit.aes = FALSE, stroke = 5) +
    geom_sf(data = dipx, inherit.aes = FALSE, fill = "pink", 
            alpha = 0.5, color = "purple3") +
    geom_point(data = morning_sleep_site, aes(x=mean_x, y=mean_y), inherit.aes = FALSE,  colour = "slateblue4", shape = 4)+
    geom_text(data = dipx_centroids, inherit.aes = FALSE,
              aes(x = X, y = Y, label = newID),
              size = 3, color = "violetred4") +
    geom_text(data = hotspot_labels,
              aes(x = mean_x, y = mean_y, label = label),
              inherit.aes = FALSE,
              color = "black", size = 7, hjust = 0, vjust = -0.4) +
    xlim(x_range) +
    ylim(y_range)+
    theme_classic() +
    labs(title = date)+
    theme(legend.position = "none",
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          title = element_text(size = 18))
  
  
  # Filter vedba data for the day
  coati_vedba_i <- coati_vedba[as.Date(coati_vedba$timestamp) == date,]
  max_vedba_i <- max(coati_vedba_i$vedba, na.rm = TRUE)
  
  # VEDBA time series plots per cluster
  cluster_ids <- sort(unique(hotspot_visits$cluster))
  cluster_plots <- list()
  
  for (c in cluster_ids) {
    hotspot_i <- hotspot_visits[hotspot_visits$cluster == c, ]
    vedba_cluster_data <- data.frame()
    
    for (j in unique(hotspot_i$ID)) {
      hotspot_c_j <- hotspot_i[hotspot_i$ID == j, ]
      vedba_j <- coati_vedba[coati_vedba$ID == j, ]
      
      vedba_i_j <- vedba_j[
        vedba_j$timestamp >= hotspot_c_j$entry_time &
          vedba_j$timestamp <= hotspot_c_j$exit_time, ]
      
      vedba_i_j$cluster <- c
      vedba_i_j$ID <- j
      vedba_cluster_data <- rbind(vedba_cluster_data, vedba_i_j)
    }
    
    p <- ggplot(vedba_cluster_data, aes(x = timestamp, y = vedba, color = ID)) +
      geom_line() +
      labs(title = paste("Cluster", c), x = "Time", y = "VEDBA") +
      ylim(0, max_vedba_i) +
      scale_color_manual(values = id_colors) +
      theme_classic() +
      geom_hline(aes(yintercept = thresh), color = "red", linetype = "dashed", linewidth = 0.8) +
      theme(legend.position = "none",  
            axis.text.x = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            title = element_text(size = 16)
            )
    
    cluster_plots[[as.character(c)]] <- p
  }
  
  # Create shared legend
  legend_plot <- ggplot(data.frame(ID = names(id_colors), vedba = 1),
                        aes(x = ID, y = vedba, color = ID)) +
    geom_point(size = 5) +
    scale_color_manual(values = id_colors) +
    theme_void() +
    theme(legend.position = "right")+
    theme(axis.text=element_text(size=30), 
          legend.text = element_text(size=25),
          legend.title = element_text(size=25))
  
  manual_legend <- cowplot::get_legend(legend_plot)
  
  # Combine VEDBA plots
  vedba_stack <- wrap_plots(cluster_plots, ncol = 1)
  final_vedba <- vedba_stack | manual_legend
  
  # Combine with map
  final_plot <- map_plot / final_vedba + plot_layout(heights = c(2, length(cluster_ids)))
  
  final_plot <- map_plot | vedba_stack | manual_legend
  
  # Optional: adjust widths if needed
  final_plot <- final_plot + plot_layout(widths = c(5, 2, 0.5))
  
  # Save
  ggsave(filename = paste0(plot_dir, 'hotspots/vedba/vedba_combined_', date, '.png'), 
         plot = final_plot, width = 30, height = 4 + 3 * length(cluster_ids), dpi = 300)
}


#for each hotspot, get the proportion of time each individual was in high vs low activity


#time of entry







