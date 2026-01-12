#this script is for looking at the trees the coatis visit when they were ripe (just the study trees where samples were taken)
#using the simulated 5 minute data and the csv sheet from fruit tree collection

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
library(fuzzyjoin)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/nutriscape_results/level1/'
gps_file <- "presidente2025_xy_10min_level0.RData"
id_file <- 'presidente2025_coati_ids.RData'
sim_5min <- 'simulated_gps_5min.RData'
ripe <- 'dipx_ripeness.RData'

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)
load(sim_5min)
load(ripe)
dipx <- st_read("C:/Users/egrout/Dropbox/coatithon/processed/2025/shapefile_data/DipteryxTreeCrowns20242025_upd1202.shp")
dipx <- st_make_valid(dipx)


#how many times each coati visits each tree:

#make gps in an sf object

dipx_utm_buffer <- st_buffer(dipx, dist = 15)

sim_sf <- st_as_sf(sim_df, coords = c("x", "y"), crs = st_crs(dipx))

#Use buffered crowns for detecting visits - this code is finding the times when the coatis overlap the dipteryx trees
sim_sf <- st_join( #st_join matches points (animal locations) with polygons (tree buffers)
  sim_sf,
  dipx_utm_buffer %>%
    select(newID) %>%
    rename(tree_id = newID),
  join = st_within #the join keeps only matches where the animal point is inside a buffered tree crown
)

#Arrange by individual & datetime
sim_sf <- sim_sf %>%
  arrange(ID, time)

# Identify when an individual enters/exits a tree
sim_sf <- sim_sf %>%
  group_by(ID) %>%
  mutate(
    visit_start = ifelse(
      !is.na(tree_id) & (is.na(lag(tree_id)) | lag(tree_id) != tree_id),
      1, 0
    ),
    event_id = ifelse(!is.na(tree_id), cumsum(visit_start), NA)   # unique visit number per individual
  ) %>%
  ungroup()

#calculate the entry and exit times for each individual into each tree
entry_exit <- sim_sf %>%
  filter(!is.na(event_id)) %>%   # keep only rows with event_id (animal in tree)
  group_by(ID, tree_id, event_id) %>%
  summarise(
    entry_time = min(time, na.rm = TRUE),
    exit_time  = max(time, na.rm = TRUE),
    duration   = as.numeric(difftime(exit_time, entry_time, units = "mins")),
    .groups = "drop"
  )

#replace duration for visit less than 5 mins to 2.5mins
entry_exit <- entry_exit %>%
  mutate(
    duration = as.numeric(difftime(exit_time, entry_time, units = "mins")),
    duration = ifelse(duration == 0, 2.5, duration)
  )

#merge ripe dataframe to entry_exit df
entry_exit$Date <- as.Date(entry_exit$entry_time)
#issue with the merge is that the same tree has multiple ripeness scores

#plotting count of ripeness levels for each tree regardless of date
p <- ggplot(ripe, aes(x = as.factor(treeID), fill = Ripeness))+
  geom_bar(stat = "count",
           position = "stack")+
  facet_wrap(~Date, scales = "free_x")+
  xlab("Tree ID")+
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 310, hjust = 0))

p
ggsave(filename = paste0(plot_dir, "ripeness_overtime.png"),
       plot = p, height = 12, width = 17)

#proportion of each tree which is unripe-riper-ripe-overripe
prop_ripe <- ripe %>%
  group_by(Date, treeID) %>%
  summarise(
    unripe = mean(Ripeness == "unripe", na.rm = TRUE),
    riper  = mean(Ripeness == "riper", na.rm = TRUE),
    ripe   = mean(Ripeness == "ripe", na.rm = TRUE),
    overripe = mean(Ripeness == "overripe", na.rm = TRUE),
    .groups = "drop"
  )

#max ripe value per fruit collection for each tree
ripe <- ripe %>%
  group_by(Date, treeID) %>%
  mutate(
    max_Ripeness = Mode(Ripeness)
    )%>%
  ungroup()

#cut duplicate rows to the unique ones
ripe <- subset(ripe, select=-c(Ripeness, fruitID))
ripe <- ripe[!duplicated(ripe), ]

#find usual time between each trees samples taken
ripe <- ripe %>%
  group_by(treeID) %>%
  arrange(Date, .by_group = TRUE) %>%  # make sure dates are in order
  mutate(
    date_diff = as.numeric(difftime(Date, lag(Date), units = "days"))
  )

g <- ggplot(ripe[!is.na(ripe$date_diff),], aes(x= as.factor(treeID), y= date_diff, group = as.factor(treeID)))+
  geom_boxplot()+
  theme_classic()+
  ylab("Duration between fruit collected (days)")+
  xlab("tree ID")+
  theme(axis.text.x = element_text(angle = 310, hjust = 0))

ggsave(filename = paste0(plot_dir, "timediffbetweensamples.png"),
       plot = g, height = 5, width = 7, dpi = 900)

#this plot shows that many trees have fruits collected once, so I think instead of a buffer, I will just match them and see if there is a general pattern with duration of visit and ripeness

#add numbers to the ripeness values (0 = unripe, 4 = overripe)
ripe$ripe_score <- NA
ripe$ripe_score[which(ripe$max_Ripeness == "unripe")] <- 1
ripe$ripe_score[which(ripe$max_Ripeness == "riper")] <- 2
ripe$ripe_score[which(ripe$max_Ripeness == "ripe")] <- 3
ripe$ripe_score[which(ripe$max_Ripeness == "overripe")] <- 4

#plot how the proportion of these categories changes over time
#using prop_ripe df as ripe df just has the max values of each ripeness for each tree, and I want to see how this proportion changes over time between the ripeness catagories
prop_ripe_long <- prop_ripe %>%
  pivot_longer(
    cols = c(unripe, riper, ripe, overripe),
    names_to = "Ripeness",
    values_to = "Proportion"
  )
#filter to when there are atleast 3 unique dates with data
trees_with_enough_data <- prop_ripe_long %>%
  group_by(treeID) %>%
  filter(n_distinct(Date) >= 3) %>%
  ungroup()

p2 <- ggplot(trees_with_enough_data, aes(x = Date, y = Proportion, color = Ripeness)) +
  geom_smooth(size = 1) +
  #facet_wrap(~treeID, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Date",
    y = "Proportion",
    title = "Ripeness Proportions Over Time"
  ) +
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 310, hjust = 0))
p2

ggsave(filename = paste0(plot_dir, "prop_ripeness_overtime.png"),
       plot = p2, height = 8, width = 10, dpi = 900)


#-----------------------------------------------------------------------
#correlate number of dipteryx visits by proportion of ripeness values
#we'd expect that there will be more visits when the proportion of ripe fruits is higher

#get mean number of dipx tree visits per ind
mean_dipx_count <- merge_df %>%
  group_by(Date, ID) %>%
  summarise(
    n_visits = n(),
    avg_duration_mins_ind = mean(duration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Date) %>%
  summarise(
    mean_n_visits = mean(n_visits, na.rm = TRUE), # mean visits per individual
    avg_duration_mins = mean(avg_duration_mins_ind, na.rm = TRUE), # mean duration per individual
    .groups = "drop"
  )

ggplot(mean_dipx_count, aes(Date, mean_n_visits))+
  geom_point()+
  ggtitle("Mean Number of Dipteryx visits per day")+
  ylab("Count")+
  geom_smooth(se = T, method = "gam", formula = y ~ s(log(x)))+
  theme_classic()

#correlate with proportion of ripe fruit for each day (regardless of tree ID)
prop_ripe_perday <-  prop_ripe %>%
  group_by(Date) %>%
  summarise(
    unripe   = mean(unripe, na.rm = TRUE),
    riper    = mean(riper, na.rm = TRUE),
    ripe     = mean(ripe, na.rm = TRUE),
    overripe = mean(overripe, na.rm = TRUE),
    .groups = "drop"
  )

ripe_dipxvisit <- merge(x = mean_dipx_count, y = prop_ripe_perday, by = "Date", all.x = T)
ripe_dipxvisit_long <- ripe_dipxvisit %>%
  pivot_longer(
    cols = c(unripe, riper, ripe, overripe),
    names_to = "Ripeness",
    values_to = "Proportion"
  )

f <- ggplot(ripe_dipxvisit_long, aes(x = mean_n_visits, y = Proportion, group = Ripeness, colour = Ripeness))+
         geom_point()+
  geom_smooth(se = T, method = "gam", formula = y ~ s(log(x)))+
  ylab("Proportion of visited Dipteryx trees")+
  xlab("Mean number of Dipteryx visits")+
  theme_classic(base_size = 15)

ggsave(filename = paste0(plot_dir, "ripeness_n_visits.png"),
       plot = f, height = 8, width = 10, dpi = 900)


#------------------------------------------------------------------

#now want to add tree ripeness for the trees that were visited that we have samples for to see if they spend more time in trees which had more ripe fruit - basically joining entry_exit df with prop_ripe df

# Ensure consistent column names for joining
ripe <- ripe %>%
  rename(tree_id = treeID)
#to increase sample size, will add the ripeness values to 1 day before and after the date the fruit sample was collected, so there is more chance of a tree visit 


ripe_extrapolate <- ripe %>%
  # make sure your tree_id column is correct
  group_by(Date, tree_id) %>%
  ungroup() %>%
  # create a column with date offsets -1, 0, +1
  tidyr::expand_grid(offset = c(-1, 0,1)) %>%
  # adjust the Date by the offset
  mutate(Date = Date + days(offset)) %>%
  select(-offset)

merge_df <- merge(x = entry_exit, y = ripe_extrapolate, by = c("Date", "tree_id"), all.x = TRUE)
merge_df_filt <- merge_df[!is.na(merge_df$max_Ripeness),] #filter to just the trees which have a ripeness value

#get the mean duration all inds visited each tree each day for plotting
merge_df_summary <- merge_df_filt %>%
  group_by(Date, tree_id, ID) %>%
  summarise(
    total_duration = sum(duration, na.rm = TRUE),
    max_Ripeness = first(max_Ripeness, default = NA),
    .groups = "drop"
  ) %>%
  group_by(Date, tree_id) %>%
  summarise(
    mean_duration = mean(total_duration, na.rm = TRUE),
    max_Ripeness = first(max_Ripeness, default = NA),  # NA if missing
    .groups = "drop"
  )


p3 <- ggplot(merge_df_summary, aes(x= Date, y= mean_duration, color = max_Ripeness))+
  geom_point(size = 3)+
  theme_classic(base_size = 15)+
  ylab("Mean Duration in each tree per day across individuals (minutes)")+
  facet_wrap(~tree_id)

ggsave(filename = paste0(plot_dir, "ripeness_durintree.png"),
       plot = p3, height = 12, width = 19, dpi = 900)



ggplot(merge_df_summary, aes(x= max_Ripeness, y= mean_duration, color = max_Ripeness))+
  geom_point(size = 3)+
  theme_classic(base_size = 15)+
  ylab("Mean Duration in each tree per day across individuals (minutes)")



#------------------------------------------------------------------
#
#which trees were the study trees? - let's plot them!

#filter the dipx to the study trees which have ripeness values
dipx_filt <- dipx[which(dipx$newID %in% unique(ripe$tree_id)),]

for(day in unique(sim_sf$date)){
  
  #day <- "2025-02-16"
  sim_df_day <- sim_df[sim_df$date == day,]
  
  # Get bounding box of days tracks
  bbox <- st_bbox(sim_sf[sim_sf$date == day,])
  
  #get the ripeness values of the study trees tracked for this day
  dipx_filt_day <- dipx %>%
    filter(newID %in% unique(ripe_extrapolate$tree_id[ripe_extrapolate$Date == day])) %>%
    left_join(
      ripe_extrapolate %>% filter(Date == day),
      by = c("newID" = "tree_id")
    )

p4 <- ggplot() +
  geom_sf(data = dipx, fill = "plum1", alpha = 0.4, color = "plum3") +
  geom_sf(data = dipx_filt, fill = "plum1", alpha = 0.7, color = "red") +
  geom_path(data = sim_df_day,
            aes(x = x, y = y, group = ID, color = ID),
            linewidth = 0.5, alpha = 0.4) +
  geom_sf(data = dipx_filt_day, aes(fill = max_Ripeness), color = "black", alpha = 0.6) +
  coord_sf(
    xlim = c(bbox$xmin, bbox$xmax),
    ylim = c(bbox$ymin, bbox$ymax)
  ) +
    scale_fill_manual(values = c(
      "unripe" = "yellow",
      "riper" = "olivedrab1",
      "ripe" = "limegreen",
      "overripe" = "green4"
    ), na.value = "grey80")+
  theme_classic(base_size = 15) +
  labs(title = paste0("Ripeness of Study Trees on ", as.Date(day)))

  
  ggsave(filename = paste0(plot_dir, "/treevisits/studytreevisits_", as.Date(day), ".png"),
         plot = p4, height = 10, width = 12, dpi = 900)
}



