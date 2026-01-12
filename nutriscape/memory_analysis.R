#This script is looking at which trees are visited and when, to determine how long the group wait until they return, aim to incorporate the fruit tree quality to see whether they just optimise visits to highly yielding trees or maximise fruit yielding success by waiting until there will be enough fruit there again.

library(cocomo)
library(dplyr)
library(sf)
library(ctmm)
library(data.table)
library(ggplot2)
library(ggrepel)  # for non-overlapping labels in ggplot
library(hms)
library(RColorBrewer)


#-------ggplot2#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level1/'
gps_file <- "presidente2025_xy_10min_level0.RData"
id_file <- 'presidente2025_coati_ids.RData'

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

dipx <- st_read("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")
dipx <- st_transform(dipx, crs = 4326)

crowns <- read.csv("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/CoatiTreeCrowns_2024-12-10_2025-02-16.csv")


#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}

#-----MAIN------
#Fist getting the UTM data cleaned and with coati names into a dataframe

breaks <- (which(diff(as_date(ts)) == 1)) + 1 #adding 1 as it gets the last each date and we need the index of the first time point for each day

clean <- preprocess_gps_level0_to_level1(xs = xs, ys = ys, timestamps = ts, ids = coati_ids, breaks = breaks, max_isolated_point_dist = 85, max_dist_percentile = 0.98, max_speed_percentile = 0.999, verbose = F)

clean_UTM_df <- matrix_to_df(xs = clean$xs, ys = clean$ys, ts = clean$timestamps)

#Get coati names in df rather than number 
#remove Ardern as not enough GPS points
coati_ids$ID <- c(1:16)
clean_UTM_df <- clean_UTM_df[!clean_UTM_df$ID == 1, ]
clean_UTM_df <- clean_UTM_df %>%
  left_join(coati_ids[, c("ID", "name")], by = "ID")
clean_UTM_df <- clean_UTM_df[as.Date(clean_UTM_df$datetime) > "2024-12-16",]
colnames(clean_UTM_df) <- c("timestamp","tag.local.identifier", "UTM_X", "UTM_Y", "location.long","location.lat","individual.local.identifier")

#made simulated gps data in first_look_at_data
load(paste0(data_dir, "simulated_gps_5min.RData"))
#Get sim_df into the same format at clean_UTM_df
sim_clean <- sim_df

sim_clean <- sim_clean %>%
  rename(UTM_X = x, UTM_Y = y) %>%
  mutate(
    timestamp = as.POSIXct(t, origin = "1970-01-01", tz = "UTC"),
    individual.local.identifier = as.character(ID),
    tag.local.identifier = as.integer(factor(individual.local.identifier))
  )

sim_sf <- st_as_sf(sim_clean, coords = c("UTM_X", "UTM_Y"), crs = 32617)
sim_sf <- st_transform(sim_sf, 4326) #get lat/lon

#Extract lon/lat back to columns
coords <- st_coordinates(sim_sf)
sim_clean$location.long <- coords[,1]
sim_clean$location.lat  <- coords[,2]

#Keep same columns and order as clean_UTM_df
sim_clean <- sim_clean %>%
  select(timestamp, tag.local.identifier, UTM_X, UTM_Y,
         location.long, location.lat, individual.local.identifier)



#---------------------------------------------------------------

#look at the time of entry and exit times into the dipteryx trees

# --- 1. Fix geometries and transform to UTM ---
dipx <- st_make_valid(dipx)
dipx_utm <- st_transform(dipx, crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# --- 2. Create a buffered version for visit detection (10 m buffer) ---
dipx_utm_buffer <- st_buffer(dipx_utm, dist = 15)

# --- 3. Convert coati data to sf ---
coati_sf <- st_as_sf(clean_UTM_df, coords = c("UTM_X", "UTM_Y"), crs = st_crs(dipx_utm))

# --- 4. Use buffered crowns for detecting visits ---
coati_sf <- st_join(
  coati_sf,
  dipx_utm_buffer %>%
    select(newID) %>%
    rename(tree_id = newID),
  join = st_within
)

# --- 5. Order by individual and time ---
coati_sf <- coati_sf %>%
  arrange(individual.local.identifier, timestamp)

# --- 6. Group consecutive points in same tree for same individual ---
dt <- as.data.table(coati_sf)
dt[, visit_group := rleid(individual.local.identifier, tree_id)]

# --- 7. Summarise to get entry/exit times ---
visits <- dt[!is.na(tree_id), .(
  entry_time = min(timestamp, na.rm = TRUE),
  exit_time  = max(timestamp, na.rm = TRUE),
  n_points   = .N
), by = .(individual.local.identifier, tree_id, visit_group)]

# --- 8. Adjust short same-block visits ---
visits <- visits %>%
  mutate(
    entry_floor = floor_date(entry_time, "5 minutes"),
    exit_floor  = floor_date(exit_time, "5 minutes"),
    exit_time   = if_else(entry_floor == exit_floor, entry_time, exit_time)
  ) %>%
  select(-entry_floor, -exit_floor, -visit_group)

# --- 9. Merge consecutive visits within 5 minutes ---
visits_dt <- as.data.table(visits)
setorder(visits_dt, individual.local.identifier, tree_id, entry_time)

visits_dt[, merge_group := cumsum(
  c(TRUE, diff(entry_time) > minutes(5) |
      tree_id[-1] != tree_id[-.N] |
      individual.local.identifier[-1] != individual.local.identifier[-.N])
)]

visits <- visits_dt[, .(
  entry_time = min(entry_time),
  exit_time  = max(exit_time),
  n_points   = sum(n_points)
), by = .(individual.local.identifier, tree_id, merge_group)]

visits[, merge_group := NULL]



#-----------------------------------------------------

plot_dir1 <- "C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level1/durintree_map/"


for (i in seq_along(unique(as.Date(coati_sf$timestamp)))) {
  
  date_i <- as.Date(unique(as.Date(coati_sf$timestamp))[i])
  
  for (ind in seq_along(unique(coati_sf$individual.local.identifier))) {
    
    individual <- unique(coati_sf$individual.local.identifier)[ind]
    
    # Filter GPS points for that day
    coati_day <- coati_sf %>%
      filter(as.Date(timestamp) == date_i,
             individual.local.identifier == individual)
    
    if (nrow(coati_day) == 0) next  # Skip if no data
    
    # Filter visits for the same date & individual
    visits_day <- visits %>%
      filter(as.Date(entry_time) == date_i,
             individual.local.identifier == individual)
    
    if (nrow(visits_day) == 0) next  # Skip if no visits
    
    # Entry & exit points
    entry_points <- coati_day %>%
      semi_join(visits_day, by = c("timestamp" = "entry_time",
                                   "individual.local.identifier"))
    exit_points <- coati_day %>%
      semi_join(visits_day, by = c("timestamp" = "exit_time",
                                   "individual.local.identifier"))
    
    # Bounding box of track
    track_bbox <- st_bbox(coati_day) %>% st_as_sfc()
    
    # All Dipteryx inside track extent
    dipx_day <- dipx_utm %>%
      filter(as.logical(st_within(geometry, track_bbox, sparse = FALSE)))
    
    # Which trees were visited
    visited_tree_ids <- unique(visits_day$tree_id)
    
    # Labels for visited trees
    dipx_day_labels <- dipx_day %>%
      left_join(
        visits_day %>%
          mutate(duration_mins = as.numeric(difftime(exit_time, entry_time, units = "mins"))) %>%
          group_by(tree_id) %>%
          summarise(duration_mins = sum(duration_mins, na.rm = TRUE), .groups = "drop") %>%
          mutate(duration_label = ifelse(duration_mins == 0, "< 5 min",
                                         paste0(round(duration_mins, 1), " min"))),
        by = c("newID" = "tree_id")
      ) %>%
      filter(!is.na(duration_mins)) %>%
      mutate(centroid = st_centroid(geometry)) %>%
      st_as_sf()
    
    # Extract coords for path
    coati_coords <- coati_day %>%
      st_coordinates() %>%
      as.data.frame() %>%
      bind_cols(coati_day %>% st_drop_geometry() %>% select(timestamp))
    
    # Plot
    gg <- ggplot() +
      geom_sf(data = dipx_day, fill = "lightgreen", alpha = 0.4, color = "darkgreen") +
      geom_sf(data = coati_day, aes(color = timestamp), size = 1) +
      geom_sf(data = entry_points, color = "black", size = 2, shape = 24, alpha = 0.5) +
      geom_path(data = coati_coords, aes(X, Y, color = timestamp), linewidth = 0.3) +
      geom_text_repel(
        data = dipx_day_labels,
        aes(
          x = st_coordinates(centroid)[, 1],
          y = st_coordinates(centroid)[, 2],
          label = paste0("Tree ", newID, "\n", duration_label)
        ),
        max.overlaps = 20,
        size = 2,
        color = "black"
      ) +
      theme_classic() +
      labs(title = paste("Coati Visits on", date_i, "-", individual)) +
      scale_color_datetime(
        date_labels = "%H:%M",
        name = "Time",
        low = "pink",
        high = "red3"
      ) +
      theme(legend.position = "right")
    
    ggsave(filename = paste0(plot_dir1, date_i, "_", individual, '.png'),
           plot = gg, width = 9, height = 9, dpi = 300)
  }
}



#----------------------------------------------------

#plot each tree on y axis against average duration of visits over time 

#duration for each visit in minutes
visits <- visits %>%
  mutate(duration_mins = as.numeric(difftime(exit_time, entry_time, units = "mins")))

visits$duration_mins[visits$duration_mins == 0] <- 2.5

visits$date <- as.Date(visits$entry_time)
# Calculate average duration per tree per day

#if visit entry or exit time is after 18:00, remove from duration calculation as they were likely a sleep tree visit

visits <- visits[format(visits$entry_time, "%H:%M:%S") < "18:00:00", ]
visits <- visits[format(visits$exit_time, "%H:%M:%S") < "18:00:00", ]

avg_tree_visits <- visits %>%
  group_by(tree_id, date) %>%
  summarise(
    avg_duration_mins = mean(duration_mins, na.rm = TRUE),
    median_duration_mins = median(duration_mins, na.rm = TRUE),
    max_duration_mins = max(duration_mins, na.rm = TRUE),
    n_visits = n(),
    .groups = "drop"
  )

#if the average visit duration is 2.5 mins, remove from plot as its likely the animal just passing through

hist(avg_tree_visits$avg_duration_mins, breaks = 100)

avg_tree_visits_filtered <- avg_tree_visits[avg_tree_visits$avg_duration_mins > 2.5, ]

n_trees <- length(unique(avg_tree_visits$tree_id))

# Make sure avg_tree_visits exists (from previous step)
gg <- ggplot(avg_tree_visits_filtered, aes(x = date, y = factor(tree_id), fill = avg_duration_mins), na.rm = T) +
  geom_tile(color = "white",  na.rm = T) +                # Each tile represents tree/day
  scale_fill_viridis_c(option = "plasma",    # Nice color scale for duration
                       name = "Mean Duration\n(minutes)") +
  labs(x = "Date", y = "Tree ID", title = "Mean Coati Visit Duration per Tree per Day") +
  geom_hline(yintercept = seq(0.5, n_trees + 0.5, by = 1), 
             color = "lightgrey", linetype = "11",  # shorter dashes and gaps
             size = 0.3) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
gg

ggsave(filename = paste0(plot_dir, 'avg_dur_dipx_pertree_with10mbuffer.png'), plot = gg, width = 9, height = 14, dpi = 300)


#----------------------------------------------------------------

visits_perday <- visits %>%
  group_by(date, individual.local.identifier) %>%
  summarise(
    n_visits = n(),
    avg_duration_mins_ind = mean(duration_mins, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(date) %>%
  summarise(
    mean_n_visits = mean(n_visits, na.rm = TRUE), # mean visits per individual
    avg_duration_mins = mean(avg_duration_mins_ind, na.rm = TRUE), # mean duration per individual
    .groups = "drop"
  )

#plot showing the mean duration of visits for all inds per day
g1 <- ggplot(visits_perday, aes(date, avg_duration_mins))+
  geom_point()+
  ylab("Mean Dipteryx Visit Duration (minutes)")+
  geom_smooth(method = "lm", se = T)+
  theme_classic()
g1

ggsave(filename = paste0(plot_dir, 'visit_dur_time_10mbuffer.png'), plot = g1, width = 7, height = 7, dpi = 300)

g2 <- ggplot(visits_perday, aes(date, mean_n_visits))+
  geom_point()+
  ggtitle("Mean Number of Dipteryx visits per day")+
  ylab("Count")+
  geom_smooth(se = T, method = "gam", formula = y ~ s(log(x)))+
  theme_classic()
g2

#?mgcv::gam

ggsave(filename = paste0(plot_dir, 'mean_visit_counts_10mbuffer.png'), plot = g2, width = 7, height = 7, dpi = 300)

hist(table(avg_tree_visits$tree_id), breaks = 30)
#cut tree visits to cases where there are at least 5 visits

avg_visits_cut <- avg_tree_visits %>%
  group_by(tree_id) %>%
  filter(n() >= 5) %>%
  ungroup()

ggplot(avg_visits_cut, aes(date, avg_duration_mins, color = as.factor(tree_id)))+
  geom_point()+
  geom_line()+
  theme_classic()+
  theme(legend.position = 'none')


#for each tree, which individuals are visiting the longest?
n_visits_each_tree <- visits %>%
  group_by(tree_id, individual.local.identifier) %>%
  summarise(
    sum_dur_visits = sum(duration_mins, na.rm = TRUE),
    .groups = "drop"
  )

hist(n_visits_each_tree$sum_dur_visits)
max(n_visits_each_tree$sum_dur_visits)

#remove Gandhi from plot as she didn't have 5/10 min resolution gps data
n_visits_each_tree <- n_visits_each_tree[!n_visits_each_tree$individual.local.identifier == "Gandhi",]
n_visits_each_tree <- n_visits_each_tree[!n_visits_each_tree$individual.local.identifier == "Einstein",] #Eisntein's GPS failed within 2 weeks

#most visited tree is 1721
ggplot(n_visits_each_tree[n_visits_each_tree$tree_id == 1721,], aes(x=as.factor(tree_id), y=sum_dur_visits, fill = individual.local.identifier))+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_discrete(name = "ID") +
  labs(
    x = "Tree ID",
    y = "Total visit duration")+
  theme_classic()

#any tree with fewer than 5 visits across all individuals is removed
tree_ind_counts <- n_visits_each_tree %>%
  filter(tree_id != 1721) %>%   # optional: remove tree 1721
  group_by(tree_id, individual.local.identifier) %>%
  summarise(total_dur_visits = sum(sum_dur_visits), .groups = "drop") %>%
  # Keep only trees with at least 20 minutes they were visited (across all individuals)
  group_by(tree_id) %>%
  filter(sum(total_dur_visits) >= 20) %>%
  ungroup()


# Make a heatmap
g3 <- ggplot(tree_ind_counts, aes(x = individual.local.identifier, y = factor(tree_id), fill = total_dur_visits)) +
  geom_tile(color = "white") +   # each cell is one tree x individual
  scale_fill_viridis_c(option = "magma", name = "Visits Duration (min)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # readable x labels
    axis.text.y = element_text(size = 5)
  ) +
  labs(
    x = "Individual ID",
    y = "Tree ID",
    title = "Total duration each coati visited each tree"
  )

ggsave(filename = paste0(plot_dir, 'visit_counts_perind.png'), plot = g3, width = 6, height = 10, dpi = 300)








