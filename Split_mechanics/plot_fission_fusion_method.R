#making plots for paper to illustrate how we extracted fission and fusion events

#LIBRARY
library(lubridate)
library(scales)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(gganimate)

#----------PARAMETERS - MODIFY HERE--------------

#which group (galaxy or presedente)
group <- 'galaxy'

#who is using (ari or emily)
user <- 'emily'

#whether to identify splits and merges automatically (if F) or use manually identified events (if T)
use_manual_events <- F

#---PARAMETERS (probably don't modify)---

#events filename - where to get the split/merge events for manually labeled events
events.filename <- paste0('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/',group,'_manual_split_merge_clean.csv') 

#radii to use
R_inner <- 15
R_outer <- 50

#------TIME ZONE------
#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS

codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
if(group == 'galaxy'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    plotdir <- "C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/"
 } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    plotdir <- "C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level2/"
 }

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(groupdir)

#read in coati ids
load(file=paste0(group,'_coati_ids.RData'))

#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})

#read in timestamp data
load(file=paste0(group,'_xy_highres_level2.RData'))

#events made with the identify_splits_and_merges_new in charecterize_splits_and_merges with level2 data
if(group == 'galaxy'){
  load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/gal_events.RData") 
} else if(group == 'presedente'){
  load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/pres_events.RData") 
}


#look at diadic distances:

#filter to shorter time window
start_day <- which(ts == "2022-01-06 11:00:00 UTC")
end_day <- which(ts == "2022-01-06 13:59:00 UTC")

xs_filt <- xs[,start_day:end_day]
ys_filt <- ys[,start_day:end_day]
ts_filt <- ts[start_day:end_day]

#put into df for easier handling
xs_df <- as.data.frame(xs_filt)
ys_df <- as.data.frame(ys_filt)

colnames(xs_df) <- ts_filt
colnames(ys_df) <- ts_filt


# Add an 'Individual' column
xs_df$Individual <- 1:nrow(xs_df)
ys_df$Individual <- 1:nrow(ys_df)

# Convert to long format
xs_long <- xs_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "X_UTM") 
ys_long <- ys_df %>%
  pivot_longer(cols = -Individual, names_to = "Time", values_to = "Y_UTM") 

# Merge the long data frames
utm_long <- xs_long %>%
  inner_join(ys_long, by = c("Individual", "Time"))

# Function to calculate dyadic distances
calculate_dyadic_distances <- function(data) {
  distances <- data %>%
    inner_join(data, by = "Time", suffix = c(".1", ".2")) %>%
    filter(Individual.1 < Individual.2) %>%
    mutate(Distance = sqrt((X_UTM.1 - X_UTM.2)^2 + (Y_UTM.1 - Y_UTM.2)^2)) %>%
    select(Time, Individual.1, Individual.2, Distance)
  
  return(distances)
}

# Calculate dyadic distances
dyadic_distances <- calculate_dyadic_distances(utm_long)




# Function to generate pairwise data frame for a specific individual
generate_distance_df <- function(individual_id) {
  individual_data <- dyadic_distances %>%
    filter(Individual.1 == individual_id | Individual.2 == individual_id) %>%
    mutate(Other_Individual = if_else(Individual.1 == individual_id, Individual.2, Individual.1)) %>%
    select(Time, Other_Individual, Distance) %>%
    mutate(Focal_Individual = individual_id)
  
  return(individual_data)
}

# Generate distance data for all individuals
distance_data_list <- lapply(1:nrow(xs), generate_distance_df)

# Combine all distance data into one data frame
combined_distance_data <- bind_rows(distance_data_list)

# Replace ID numbers with actual names
combined_distance_data$Focal_Individual <- coati_ids$name[combined_distance_data$Focal_Individual]
combined_distance_data$Other_Individual <- coati_ids$name[combined_distance_data$Other_Individual]

# Fill missing data using linear interpolation and carry forward/backward for leading/trailing NAs
combined_distance_data <- combined_distance_data %>%
  group_by(Focal_Individual, Other_Individual) %>%
  mutate(Distance = zoo::na.locf(zoo::na.locf(zoo::na.approx(Distance, na.rm = FALSE), na.rm = FALSE), fromLast = TRUE)) %>%
  ungroup()


com_sat <- combined_distance_data[combined_distance_data$Other_Individual == "Venus" & combined_distance_data$Focal_Individual == "Planeta",]

com_sat$UTCTime <- as.POSIXct(com_sat$Time, format = "%Y-%m-%d %H:%M:%S")
com_sat$Time <- com_sat$UTCTime  - 5 * 60 * 60 

com_sat <- com_sat %>%
  mutate(under_50 = Distance < 50,
         under_15 = Distance < 15,
         group_under_50 = cumsum(c(0, diff(as.integer(under_50))) != 0),
         group_under_15 = cumsum(c(0, diff(as.integer(under_15))) != 0))

intervals <- com_sat %>%
  filter(under_50) %>%
  group_by(group_under_50) %>%
  summarise(start_time = min(Time), end_time = max(Time))

intervals_under_15 <- com_sat %>%
  filter(under_15) %>%
  group_by(group_under_15) %>%
  summarise(start_time = min(Time), end_time = max(Time))





# Create the plot with shaded areas
g <- ggplot(com_sat, aes(x = Time, y = Distance)) +
  geom_line() +
  geom_rect(data = intervals, aes(xmin = start_time, xmax = end_time, ymin = -Inf, ymax = Inf), fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  geom_rect(data = intervals_under_15, aes(xmin = start_time, xmax = end_time, ymin = -Inf, ymax = Inf), fill = "darkblue", alpha = 0.3, inherit.aes = FALSE) +
  scale_x_datetime(date_breaks = "30 min", labels = date_format("%H:%M")) + # Adjust the datetime scale 
  geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(x = "Time",
    y = "Dyadic distance (m)",
    title = "Identifying periods of connectivity",
    subtitle = "Light shaded area indicates when distance < 50 m
Dark shaded area indicates when distance < 15 m"
  )+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16))

g

ggsave(g, file = paste0(plotdir, "diadic_dist.png"), width = 10, height = 6)




#now plotting gps points to show chain rule
for (i in seq(1, length(utm_long$Time), by = 120)) {
  
  # Filter the data for the specific time point
  utm_long_filt <- utm_long[utm_long$Time == utm_long$Time[i],]
  
  # Count number of individuals with valid X_UTM and Y_UTM data
  num_valid_individuals <- sum(!is.na(utm_long_filt$X_UTM) & !is.na(utm_long_filt$Y_UTM))
  
  # Skip to next iteration if fewer than 8 individuals have valid data
  if (num_valid_individuals < 8) next
  
  # Create a dataframe to store pairwise distances
  distances <- expand.grid(Individual1 = utm_long_filt$Individual, Individual2 = utm_long_filt$Individual) %>%
    filter(Individual1 != Individual2) %>%
    left_join(utm_long_filt, by = c("Individual1" = "Individual")) %>%
    rename(X_UTM1 = X_UTM, Y_UTM1 = Y_UTM) %>%
    left_join(utm_long_filt, by = c("Individual2" = "Individual")) %>%
    rename(X_UTM2 = X_UTM, Y_UTM2 = Y_UTM)
  
  # Calculate the Euclidean distance
  distances <- distances %>%
    mutate(distance = sqrt((X_UTM1 - X_UTM2)^2 + (Y_UTM1 - Y_UTM2)^2))
  
  # Filter pairs with distance less than 50 meters
  close_pairs <- distances %>%
    filter(distance < 50)
  
  # Calculate the extent of the plot
  x_range <- range(utm_long_filt$X_UTM, na.rm = TRUE)
  y_range <- range(utm_long_filt$Y_UTM, na.rm = TRUE)
  
  # Calculate coordinates for scale bar
  scale_bar_x <- x_range[1] + (x_range[2] - x_range[1]) * -0.1  # Adjust position along x-axis
  scale_bar_y <- y_range[1] + (y_range[2] - y_range[1]) * 0.1  # Adjust position along y-axis
  scale_bar_length <- 50  # Length of the scale bar in meters
  scale_bar_label <- "50 m"  # Label for the scale bar
  
  # Plot
  g2 <- ggplot(utm_long_filt, aes(x = X_UTM, y = Y_UTM)) +
    geom_point(size = 5, color = "purple4") +
    #geom_circle(data = utm_long_filt, aes(x0 = X_UTM, y0 = Y_UTM, r = 25), fill = NA, color = "skyblue2") + # Add circles with radius of 25 meters (diameter of 50 meters)
    geom_segment(data = close_pairs, aes(x = X_UTM1, y = Y_UTM1, xend = X_UTM2, yend = Y_UTM2), color = "red") +
    geom_point(size = 5, color = "purple4") +
    geom_segment(x = scale_bar_x, xend = scale_bar_x + scale_bar_length, y = scale_bar_y, yend = scale_bar_y, color = "black", size = 1) +
    geom_text(x = scale_bar_x + scale_bar_length / 2, y = scale_bar_y - 0.05, label = scale_bar_label, vjust = 1.5, hjust = 0.5, size = 7) +
    coord_fixed() +
    theme_classic() +
    labs(
      x = "X UTM",
      y = "Y UTM") +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = 'white', colour = 'white'),
          plot.background = element_rect(fill="white"))
  
  # Save the plot
  ggsave(g2, file = paste0(plotdir, "chainrule_plots/nocircles/", i, ".png"), width = 6, height = 6)
  
}




