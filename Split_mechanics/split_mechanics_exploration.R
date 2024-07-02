#LIBRARY
library(lubridate)
library(scales)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggtext)
library(glue)
library(tidyr)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS

#who is using (ari or emily)
user <- 'emily'

#which group - galaxy or presedente
group <- 'presedente' #subdirectory where the group data is stored

#whether to identify splits and merges automatically (if F) or use manually identified events (if T)
use_manual_events <- F

#choose whether want the male events or non-male events
with_males <- F

if(user %in% c('Ari','ari')){
  codedir <- '~/Dropbox/code_ari/coatithon_code/'
  dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
  if(group == 'galaxy'){
    groupdir <- '~/Dropbox/coati/processed/galaxy/'
  } else if(group=='presedente'){
    groupdir <- '~/Dropbox/coati/processed/presedente/'
  }
} else{
  codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
  if(group == 'galaxy'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level2/'
  }
}

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(codedir)


if(use_manual_events){
  events <- read.csv(paste0('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/',group,'_manual_split_merge_clean2.csv'),  sep=",", header=TRUE)
  #preprocess events to...
  events <- events[which(events$fission_time!='before start'),] #remove events where we missed the start
  events <- events[which(events$event_type %in% c('fission','fusion')),] #only include fission and fusion events (remove 'almost fusion')
  
  
} else{ #otherwise load in automated events
  #read in automated events - df made in characterize_splits_amd_merges 
  load(paste0(groupdir, group,'_auto_ff_events_characterized.RData'))
  
}

#read in coati ids
setwd(groupdir)

#read in timestamp data
load(file=paste0(group,'_xy_highres_level2.RData'))

load(file=paste0(group,'_coati_ids.RData'))
#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})


if(group == 'presedente'){
  high_col <- "turquoise4"
  comp_col <- "paleturquoise2"
  break_1 <- 12
  break_2 <- 40
  break_3 <- 5
  break_4 <- 20
  break_5 <- 20
  break_6 <- 50
} else if(group == "galaxy"){
  high_col <- "darkolivegreen4"
  comp_col <- "darkolivegreen1"
  break_1 <- 10
  break_2 <- 40
  break_3 <- 5
  break_4 <- 20
  break_5 <- 50
  break_6 <- 50
}

#redo matrix with and without males to see how this affects the patterns observed

#find the events where males are involved 
if(group == 'presedente'){
  male_coatis <- c("Sam", "Ken", "Gen", "Lul", "Man")
} else if(group == "galaxy"){
  male_coatis <- "Gus"
}

# Function to check if all individuals in a group are males
all_males <- function(group, male_names) {
  individuals <- strsplit(group, ",//s*")[[1]]
  all(individuals %in% male_names)
}

# Create a new dataframe based on the choice
if (with_males) {
  events <-  events %>%
    rowwise() %>%
    filter(all_males(group_A, male_coatis) | all_males(group_B, male_coatis))
  name <- "with males"
} else {
  #make event dataframe without males
  events <- events %>%
    rowwise() %>%
    filter(!(all_males(group_A, male_coatis) | all_males(group_B, male_coatis)))
  name <- "without males"
} 


#making speeds vs group size dataframe for gro_size_dist_travelled_plot
create_group_distance_df <- function(event_type) {
  # Subset dataframe based on event_type and columns
  event_df <- na.omit(events[events$event_type == event_type, c("n_A", "n_B", "A_during_disp", "B_during_disp")])
  
  # Create Distance_Larger_Group and Distance_Smaller_Group columns
  event_df$Distance_Larger_Group <- ifelse(event_df$n_A >= event_df$n_B, event_df$A_during_disp, event_df$B_during_disp)
  event_df$Distance_Smaller_Group <- ifelse(event_df$n_A < event_df$n_B, event_df$A_during_disp, event_df$B_during_disp)
  
  save(event_df, file = paste0(groupdir, group, "_", event_type, "_grpsize.RData"))
  
  # Return the dataframe
  return(event_df)
}

# Usage for "fusion" event type
fusion_df <- create_group_distance_df("fusion")
# Usage for "fission" event type
fission_df <- create_group_distance_df("fission")



#Deciding cut-off points to categorise different types of fissions and fusions

# Combine the two columns into a single vector, removing NA values
during_dist <- as.vector(t(rbind(events$A_during_disp, events$B_during_disp)))
during_dist <- na.omit(during_dist)

# Compute density
density_data <- density(during_dist)

# Create a histogram with density on the y-axis
hist(during_dist, breaks = 100, freq = FALSE, main = " ", xlab = "Distance during displacement", ylab = "Density", col = "darkslategray4")
lines(density(during_dist), col = "darkorange", lwd = 2)

#Finding the point where to split the data to slow group and fast group (by finding the lowest peak between the high peaks)
yvals <- density(during_dist)$y
xvals <- density(during_dist)$x
d_yvals <- yvals[2:length(yvals)] - yvals[1:length(yvals)-1]
d_xvals <- xvals[2:length(xvals)] - xvals[1:length(xvals)-1]
d <- d_yvals/d_xvals

d <- d[(15 < xvals) & (xvals < 30)]
xvals <- xvals[(15 < xvals) & (xvals < 30)]

index_min <- which(diff(sign(d)) != 0)
threshold <- xvals[index_min]
print(threshold)
abline(v=threshold, lty = 3)


# Combine the two columns into a single vector, removing NA values
diff_dist <- abs(events$A_during_disp - events$B_during_disp)
events$dist_diff <-  diff_dist

# Create a histogram with density on the y-axis
hist(diff_dist, breaks = 100, freq = FALSE, main = " ", xlab = "Distance difference during displacement", ylab = "Density", col = "darkslategray3")
lines(density(na.omit(diff_dist)), col = "darkorange", lwd = 2)
abline(v=threshold, lty = 3)

#get speed to plot the histogram for speed during displacement - m/s
events$B_during_speed <- events$B_during_disp/(events$end_time - events$start_time)
events$A_during_speed <- events$A_during_disp/(events$end_time - events$start_time)
diff_speed <- abs(events$A_during_speed - events$B_during_speed)

hist(diff_speed, breaks = 100, freq = FALSE, main = " ", xlab = "Speed difference during displacement", ylab = "Density", col = "darkslategray3")
lines(density(na.omit(diff_speed)), col = "darkorange", lwd = 2)

#distribution of speed for fissions
fiss_speed <- rbind(events$B_during_speed[events$event_type == "fission"], events$A_during_speed[events$event_type == "fission"])

#what was the speed of the full group before a fission?
events$AB_before_speed <- events$AB_before_disp/(events$start_time - events$before_time)

#get the max value for the histograms to not cut info
max <- max(cbind(max(events$AB_before_speed, na.rm = T), max(fiss_speed, na.rm = T)))

png(height = 800, width = 1200, units = 'px', filename = paste0(plot_dir,group,'_speed_before_during_fission.png'))
par(mar = c(5, 5, 2, 2)) #bottom, left, top, right
hist(events$AB_before_speed[events$event_type == "fission"], breaks = break_1, freq = FALSE, main = " ", xlab = "Speeds before and during fissions (m/s)", ylab = "Density", cex.lab = 2.5, cex.axis = 2.5, col = alpha(high_col, 1), xlim = c(0, 2.5), yaxt = "n", xaxt = "n")
axis(2, at = c(0,2,4,6,8,10,12), cex.axis = 2.5, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,0.5, 1, 1.5, 2, 2.5), cex.axis = 2.5, las = 1, pos = 0, padj = 0.8)
hist(fiss_speed, breaks = break_2, freq = FALSE, main = " ", col = alpha(comp_col,0.5), add = T)
legend("topright", legend=c("Before fission", "During fission"),
       fill=c(high_col, comp_col), lty=0, cex=2)
dev.off()


#distribution of speed after fissions
events$B_after_speed <- events$B_after_disp/(events$after_time - events$end_time)
events$A_after_speed <- events$A_after_disp/(events$after_time - events$end_time)
fiss_after_speed <- rbind(events$B_after_speed[events$event_type == "fission"], events$A_after_speed[events$event_type == "fission"])

png(height = 800, width = 1200, units = 'px', filename = paste0(plot_dir,group,'_speed_after_fission.png'))
par(mar = c(5, 5, 2, 2)) #bottom, left, top, right
hist(fiss_speed, breaks = break_2, freq = FALSE, main = " ", xlab = "Speeds after fissions (m/s)", ylab = "Density", cex.lab = 2.5, cex.axis = 2.5, col = alpha(high_col, 1), xlim = c(0, 2.5), ylim = c(0, 7), yaxt = "n", xaxt = "n")
axis(2, at = c(0,2,4,6,8,10,12), cex.axis = 2.5, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,0.5, 1, 1.5, 2, 2.5), cex.axis = 2.5, las = 1, pos = 0, padj = 0.8)
hist(fiss_after_speed, breaks = break_5, freq = FALSE, main = " ", col = alpha(comp_col,0.5), add = T)
legend("topright", legend=c("During fission", "After fission"),
       fill=c(high_col, comp_col), lty=0, cex=2)
dev.off()

#------------------------------------------------------
#having the speed plots before, during, and after in a panel instead
#need to restructure the data so I have a long dataframe with the before, during, and after speeds in the same dataframe
fiss_before_speed <- events$AB_before_speed[events$event_type == "fission"]
fiss_before_speed <- as.data.frame(t(fiss_before_speed))
fiss_before_speed$time <- "before"
fiss_speed <- rbind(events$B_during_speed[events$event_type == "fission"], events$A_during_speed[events$event_type == "fission"])
fiss_speed <- as.data.frame(fiss_speed)
fiss_speed$time <- "during"
fiss_after_speed <- rbind(events$B_after_speed[events$event_type == "fission"], events$A_after_speed[events$event_type == "fission"])
fiss_after_speed <-  as.data.frame(fiss_after_speed)
fiss_after_speed$time <- "after"
speeds_across_fission <- as.data.frame(rbind(fiss_before_speed, fiss_speed, fiss_after_speed))
#pivot longer for plotting
speeds_across_fission <- speeds_across_fission %>%
  pivot_longer(!time, names_to = "event_time", values_to = "count")
speeds_across_fission <- na.omit(speeds_across_fission)

speeds_across_fission$time <- factor(speeds_across_fission$time, levels = c("before", "during", "after"))

# Create the boxplot
ggplot(speeds_across_fission, aes(x = time, y = count, fill = high_col)) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Period", y = "Speed (m/s)", title = "Fission")+
  guides(fill="none")

#-------------------------------------------------------------------

#distribution of speed for fusions
events$B_before_speed <- events$B_before_disp/(events$start_time - events$before_time)
events$A_before_speed <- events$A_before_disp/(events$start_time - events$before_time)
fus_before_speed <- rbind(events$B_before_speed[events$event_type == "fusion"], events$A_before_speed[events$event_type == "fusion"])
fus_during_speed <- rbind(events$B_during_speed[events$event_type == "fusion"], events$A_during_speed[events$event_type == "fusion"]) 

png(height = 800, width = 1200, units = 'px', filename = paste0(plot_dir,group,'_speed_before_during_fusion.png'))
par(mar = c(5, 5, 2, 2)) #bottom, left, top, right
hist(fus_before_speed, breaks = break_4, freq = FALSE, main = " ", xlab = "Speeds before and during fusions (m/s)", ylab = "Density", cex.lab = 2.5, cex.axis = 2.5, col = alpha(high_col, 1), xlim = c(0, 2),ylim = c(0,7),  yaxt = "n", xaxt = "n")
axis(2, at = c(0,2,4,6,8,10,12), cex.axis = 2.5, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,0.5, 1, 1.5, 2, 2.5), cex.axis = 2.5, las = 1, pos = 0, padj = 0.8)
hist(fus_during_speed, breaks = break_6, freq = FALSE, main = " ", col = alpha(comp_col,0.5), add = T)
legend("topright", legend=c("Before fusion", "During fusion"),
       fill=c(high_col, comp_col), lty=0, cex=2)
dev.off()

#what was the speed of the full group after a fusion?
events$AB_after_speed <- events$AB_after_disp/(events$after_time - events$end_time)

#for galaxy - breaks are 50 and 5 
png(height = 800, width = 1200, units = 'px', filename = paste0(plot_dir,group,'_speed_after_fusion.png'))
par(mar = c(5, 5, 2, 2)) #bottom, left, top, right
hist(fus_during_speed, breaks = break_6, freq = FALSE, main = " ", xlab = "Speeds after fusions (m/s)", ylab = "Density", cex.lab = 2.5, cex.axis = 2.5, col = alpha(high_col, 1), xlim = c(0, 2),ylim = c(0,10),  yaxt = "n", xaxt = "n")
axis(2, at = c(0,2,4,6,8,10,12), cex.axis = 2.5, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,0.5, 1, 1.5, 2, 2.5), cex.axis = 2.5, las = 1, pos = 0, padj = 0.8)
hist(events$AB_after_speed[events$event_type == "fusion"], breaks = break_3, freq = FALSE, main = " ", col = alpha(comp_col,0.5), add = T)
legend("topright", legend=c("During fusion", "After fusion"),
       fill=c(high_col, comp_col), lty=0, cex=2)
dev.off()

#-----------------------------------------------------------
#having the speed plots before, during, and after in a panel instead
#need to restructure the data so I have a long dataframe with the before, during, and after speeds in the same dataframe
fus_before_speed <- rbind(events$B_before_speed[events$event_type == "fusion"], events$A_before_speed[events$event_type == "fusion"])
fus_before_speed <- as.data.frame(fus_before_speed)
fus_before_speed$time <- "before"
fus_during_speed <- rbind(events$B_during_speed[events$event_type == "fusion"], events$A_during_speed[events$event_type == "fusion"]) 
fus_during_speed <- as.data.frame(fus_during_speed)
fus_during_speed$time <- "during"
fus_after_speed <- as.data.frame(t(events$AB_after_speed[events$event_type == "fusion"]))
fus_after_speed$time <- "after"

speeds_across_fusion <- as.data.frame(rbind(fus_before_speed, fus_during_speed, fus_after_speed))
#pivot longer for plotting
speeds_across_fusion <- speeds_across_fusion %>%
  pivot_longer(!time, names_to = "event_time", values_to = "count")
speeds_across_fusion <- na.omit(speeds_across_fusion)

speeds_across_fusion$time <- factor(speeds_across_fusion$time, levels = c("before", "during", "after"))

# Create the boxplot
ggplot(speeds_across_fusion, aes(x = time, y = count, fill = comp_col)) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Period", y = "Speed (m/s)", title = "Fusion")+
  guides(fill="none")


#combine fission and fusion plot:
speeds_across_fission$event_type <- "fission"
speeds_across_fusion$event_type <- "fusion"

both_speeds <- rbind(speeds_across_fission, speeds_across_fusion)

# Create the boxplot
ggplot(both_speeds, aes(x = time, y = count)) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Period", y = "Speed (m/s)", title = group)+
  guides(fill="none")+
  facet_wrap(~event_type)


ggsave(paste0(plot_dir, "speeds_fis_fus.png"), width = 10, height = 5)

if(group == "galaxy"){
  both_speeds_galaxy <- both_speeds
  save(both_speeds_galaxy, file= paste0(groupdir, "speeds_galaxy.RData"))
} else if(group == "presedente"){
  both_speeds_pres <- both_speeds
  save(both_speeds_pres, file= paste0(groupdir, "speeds_presidente.RData"))
  
}




#decided to cut-off the "non-moving" from the "moving" group at 10 m based on visual inspection of the plots and from biological reasoning 
#Create column 'subgroup_move' using ifelse

dist_thresh <- 10
type <- "fission"

events$subgroup_move <- ifelse(events$A_during_disp > dist_thresh & events$B_during_disp > dist_thresh, "Both moved", ifelse(events$A_during_disp > dist_thresh, "Either moved",ifelse(events$B_during_disp > dist_thresh, "Either moved",  "Neither moved")))

# Create the new column to see if the distances moved are similar 
events$speed_comparison <- ifelse(abs(events$A_during_disp - events$B_during_disp) < dist_thresh, "Similar distance travelled", "Different distance travelled")

#was it a group or an individual who moved
#not sure this makes sense as it says singleton if just one of the groups contained a singleton
#events$individual_movement <- ifelse((events$A_during_disp > dist_thresh & events$n_A == 1) | 
#                                       (events$B_during_disp > dist_thresh & events$n_B == 1), 
#                                     "singleton", "multiple individuals")




during_dist_fission <- na.omit(as.vector(t(rbind(events$A_during_disp[events$event_type == "fission"], events$B_during_disp[events$event_type == "fission"]))))
# Compute density
density_data <- density(during_dist_fission)

# plot 1: Create a histogram with density on the y-axis
p1 <- ggplot() +
  geom_histogram(aes(x = during_dist_fission, y = after_stat(density)), breaks = seq(min(density_data$x), max(density_data$x), length.out = 100), fill = high_col, color = "black", alpha = 0.5) +
  geom_line(data = data.frame(x = density_data$x, y = density_data$y), aes(x = x, y = y), color = "darkorange", linewidth = 2) +
  labs(title = "", x = "Distance during fission", y = "Density")+
  geom_vline(xintercept = 10, linetype = "dashed", color = "black")+
  xlim(0, max(during_dist_fission))+
  theme_classic()
p1

file_path <- file.path(plot_dir, paste0(group, "_distance_displacement.png"))
ggsave(file_path, p1, width = 8, height = 4)



#plot 2: both groups moved
both_moved <- events[events$subgroup_move == "Both moved",]

p2 <- ggplot(na.omit(both_moved[both_moved$event_type == type,]), aes(x= speed_comparison))+
  geom_bar(position = "dodge", fill =c("yellow3", "darkslategray4") )+
  labs(title = "Both subgroups moved",
       x = " ",
       y = "Count",
       fill = " ") +
  ylim(0,22)+
  theme_minimal()

p2

#plot 3: either group moved
either_moved <- events[events$subgroup_move == "Either moved",]

p3 <- ggplot(na.omit(either_moved[either_moved$event_type == type,]), aes(x= speed_comparison))+
  geom_bar(position = "stack", fill = c("yellow3", "darkslategray4"))+
  labs(title = "Either subgroups moved",
       x = " ",
       y = "Count",
       fill = " ") +
  guides(fill="none")+
  theme_minimal()
p3


#plot 4: distribution of subgroup sizes

sub_size <- as.vector(t(rbind(events$n_A[events$event_type == "fission"], events$n_B[events$event_type == "fission"])))

p4 <- ggplot() +
  geom_histogram(aes(x = sub_size, y = after_stat(density)),bins = 16, fill = "darkslategray4", color = "white")+
  labs(title = "Number of individuals in subgroups during fission events",
       x = "Number of Individuals",
       y = "Frequency")

g <- (p1 + p4)/( p2 + p3)
g

#could split the plot 3 into the bigger vs smaller subgroup to move

#-------------------------------------------------------------------------------

#making dataframe where we classify type of events and save in a dataframe called filt_events

#Add split_type and subgroup_comp columns to the events dataframe

events <- events %>%
  mutate(
    split_type = case_when(
      event_type == "fission" & AB_before_disp < 10 & B_during_disp > 10 & A_during_disp > 10 ~ "bothstill_bothmove",
      event_type == "fission" & AB_before_disp < 10 & B_during_disp > 10 & A_during_disp < 10 ~ "bothstill_onemove",
      event_type == "fission" & AB_before_disp < 10 & B_during_disp < 10 & A_during_disp > 10 ~ "bothstill_onemove",
      event_type == "fission" & AB_before_disp > 10 & B_during_disp < 10 & A_during_disp > 10 ~ "bothmove_onemove",
      event_type == "fission" & AB_before_disp > 10 & B_during_disp > 10 & A_during_disp < 10 ~ "bothmove_onemove",
      event_type == "fission" & AB_before_disp > 10 & B_during_disp > 10 & A_during_disp > 10 ~ "bothmove_bothmove",
      event_type == "fission" & AB_before_disp < 10 & B_during_disp < 10 & A_during_disp < 10 ~ "bothstill_bothstill",
      
      event_type == "fusion" & AB_after_disp < 10 & B_during_disp > 10 & A_during_disp > 10 ~ "bothmove_bothstill",
      event_type == "fusion" & AB_after_disp < 10 & B_during_disp > 10 & A_during_disp < 10 ~ "onemove_bothstill",
      event_type == "fusion" & AB_after_disp < 10 & B_during_disp < 10 & A_during_disp > 10 ~ "onemove_bothstill",
      event_type == "fusion" & AB_after_disp > 10 & B_during_disp < 10 & A_during_disp > 10 ~ "onemove_bothmove",
      event_type == "fusion" & AB_after_disp > 10 & B_during_disp > 10 & A_during_disp < 10 ~ "onemove_bothmove",
      event_type == "fusion" & AB_after_disp > 10 & B_during_disp > 10 & A_during_disp > 10 ~ "bothmove_bothmove",
      event_type == "fusion" & AB_after_disp < 10 & B_during_disp < 10 & A_during_disp < 10 ~ "bothstill_bothstill",
      TRUE ~ "other"
    ),
    subgroup_comp = case_when(
      n_A == 1 & n_B == 1 ~ "one/one",
      n_A > 1 & n_B == 1 ~ "one/many",
      n_A == 1 & n_B > 1 ~ "one/many",
      n_A > 1 & n_B > 1 ~ "many/many",
      TRUE ~ "other"
    ),
    subgroup_moved = case_when(
      A_during_disp > 10 & B_during_disp <= 10 ~ "A",
      B_during_disp > 10 & A_during_disp <= 10 ~ "B",
      A_during_disp > 10 & B_during_disp > 10 ~ "both",
      A_during_disp <= 10 & B_during_disp <= 10 ~ "neither",
      TRUE ~ "unknown"
    )
  ) %>%
  na.omit()  # Remove rows with NA


detailed_events <- events
#save events 
if(group == 'galaxy'){
  save(detailed_events, file = paste0(groupdir, group, name, "_detailed_events.RData"))
}else if(group == "presedente"){
  save(detailed_events, file = paste0(groupdir, group,name, "_detailed_events.RData"))
}


#-------------------------------------------------------
#-------------------------------------------------------

#define type of event to look at
event_type <- "fusion"

#------------------------------------------------------
#------------------------------------------------------


event_type_df <- events[events$event_type == event_type,]

#remove rows where both are still to both are still as not really an event
event_type_df <- event_type_df[!event_type_df$split_type == "bothstill_bothstill",]
event_type_df <- event_type_df[!event_type_df$split_type == "other",]

contingency_table <- table(event_type_df$split_type, event_type_df$subgroup_comp)
contingency_matrix <- as.matrix(contingency_table)
print(contingency_matrix)

# Convert the matrix to a data frame for ggplot2
contingency_df <- as.data.frame(as.table(contingency_matrix))

# Plot the heatmap
all <- ggplot(contingency_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = high_col) +
  labs(title = paste(group, event_type, "heatmap"),
       x = "Subgroup Composition",
       y = "Event Type",
       fill = "Count") +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 0, vjust = 5, hjust = 0.5),
        axis.text.y = ggtext::element_markdown(), # Adjust margin here
        panel.grid.major = element_blank(),  # Remove major grid lines if desired
        panel.grid.minor = element_blank(),  # Remove minor grid lines if desired
        panel.border = element_blank(),      # Remove the outline around the plot
        panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
        plot.background = element_rect(fill = "white", color = NA)   # Ensure plot background is blank
  )+
  scale_y_discrete(labels = function(x) glue::glue(" <img src = 'C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/arrows/{x}.png' height = 50 /> "))
all

file_path <- file.path(plot_dir, paste(event_type, "matrix.png"))
ggsave(file_path, all, width = 10, height = 8)

#now want to look at the individual who moves to see if this is the singleton or the group 
#for events where both are still and one moved, was it the single individual who moved?

#change coat_ids sub-adults to subadults
coati_ids$age <- gsub("Sub-adult", "Subadult", coati_ids$age)

age_sex_df <- event_type_df[,c("event_idx", "event_type", "n_A", "n_B", "split_type", "subgroup_comp", "subgroup_moved")]


age_sex_df$A_Adult_Female <- NA
age_sex_df$A_Adult_Male <- NA
age_sex_df$A_Subadult_Female <- NA
age_sex_df$A_Subadult_Male <- NA
age_sex_df$A_Juvenile_Female <- NA
age_sex_df$A_Juvenile_Male <- NA
age_sex_df$B_Adult_Female <- NA
age_sex_df$B_Adult_Male <- NA
age_sex_df$B_Subadult_Female <- NA
age_sex_df$B_Subadult_Male <- NA
age_sex_df$B_Juvenile_Female <- NA
age_sex_df$B_Juvenile_Male <- NA

for (i in 1:nrow(event_type_df)){
  
  A_age_sex <- coati_ids[detailed_events$group_A_idxs[i][[1]],c("age", "sex")]
  A_age_sex <- paste("A", A_age_sex$age, A_age_sex$sex, sep = "_")
  age_sex_df[which(age_sex_df$event_idx == detailed_events$event_idx[i]),which(names(age_sex_df)%in% A_age_sex)]<- 1
  
  B_age_sex <- coati_ids[detailed_events$group_B_idxs[i][[1]],c("age", "sex")]
  B_age_sex <- paste("B", B_age_sex$age, B_age_sex$sex, sep = "_")
  age_sex_df[which(age_sex_df$event_idx == detailed_events$event_idx[i]),which(names(age_sex_df)%in% B_age_sex)]<- 1
  
}


#---------------------------------------------
#----AGE SEX CLASS PLOTS-------
# Summarize counts of each age/sex class based on subgroup_moved and event_type
age_sex_summary <- age_sex_df %>%
  group_by(event_type, split_type, subgroup_comp) %>%
  summarise(
    Adult_Female = sum(if_else(subgroup_moved %in% c("A", "both"), A_Adult_Female, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Adult_Female, 0, missing = 0), na.rm = TRUE),
    Adult_Male = sum(if_else(subgroup_moved %in% c("A", "both"), A_Adult_Male, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Adult_Male, 0, missing = 0), na.rm = TRUE),
    Subadult_Female = sum(if_else(subgroup_moved %in% c("A", "both"), A_Subadult_Female, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Subadult_Female, 0, missing = 0), na.rm = TRUE),
    Subadult_Male = sum(if_else(subgroup_moved %in% c("A", "both"), A_Subadult_Male, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Subadult_Male, 0, missing = 0), na.rm = TRUE),
    Juvenile_Female = sum(if_else(subgroup_moved %in% c("A", "both"), A_Juvenile_Female, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Juvenile_Female, 0, missing = 0), na.rm = TRUE),
    Juvenile_Male = sum(if_else(subgroup_moved %in% c("A", "both"), A_Juvenile_Male, 0, missing = 0), na.rm = TRUE) + 
      sum(if_else(subgroup_moved %in% c("B", "both"), B_Juvenile_Male, 0, missing = 0), na.rm = TRUE)
  )

# Convert to long format
age_sex_long <- age_sex_summary %>%
  pivot_longer(
    cols = starts_with("Adult_") | starts_with("Subadult_") | starts_with("Juvenile_"),
    names_to = "age_sex_class",
    values_to = "count"
  )

# PLOT
age_sex_long$age_sex_class <- sub("_", " ", age_sex_long$age_sex_class)
age_sex_long$age_sex_class <- as.factor(age_sex_long$age_sex_class)
age_sex_long$age_sex_class <- factor(age_sex_long$age_sex_class, levels = c("Adult Female", "Adult Male", "Subadult Female", "Subadult Male", "Juvenile Female", "Juvenile Male"))
age_sex_long$subgroup_comp <- factor(age_sex_long$subgroup_comp, levels = c("one/one", "one/many", "many/many"))

custom_colors <- c(
  "Adult Male" = "cadetblue3",
  "Subadult Male" = "coral3",
  "Juvenile Male" = "darkolivegreen3",
  "Adult Female" = "lightblue1",
  "Subadult Female" = "salmon",
  "Juvenile Female" = "darkolivegreen1"
)

g1 <- ggplot(age_sex_long, aes(x = split_type, y = count, fill = age_sex_class)) +
  geom_bar(stat = "identity", color = NA) +  # Ensure color = NA here
  scale_fill_manual(values = custom_colors) +  # Custom color scale if needed
  xlab(paste0(event_type, " type")) + ylab("Count") +
  scale_x_discrete(labels = function(x) glue::glue("<img src='C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/arrows/{x}.png' height='30' />")) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = ggtext::element_markdown(margin = margin(t = -10)),  # Adjust margin here
    panel.grid.major = element_blank(),  # Remove major grid lines if desired
    panel.grid.minor = element_blank(),  # Remove minor grid lines if desired
    panel.border = element_blank(),      # Remove the outline around the plot
    panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
    plot.background = element_rect(fill = "white", color = NA)   # Ensure plot background is blank
  ) +
  guides(fill = guide_legend(title = "Age/Sex class")) +
  facet_wrap(~ subgroup_comp)
g1

ggsave(g1 <- paste0(plot_dir, group,"_", event_type, "_agesex_nomales.png"), width = 11, height = 8)


#now looking at the age/sex classes for the "one" or "many" groups:

# Filter the data for "one/many" subgroup_comp
one_many_data <- age_sex_df %>%
  filter(subgroup_comp == "one/many")

# Prepare summaries for the "one" and "many" groups
one_many_counts <- one_many_data %>%
  mutate(
    group = case_when(
      subgroup_moved == "A" & n_A == 1 ~ "one",
      subgroup_moved == "A" & n_A > 1 ~ "many",
      subgroup_moved == "B" & n_B == 1 ~ "one",
      subgroup_moved == "B" & n_B > 1 ~ "many",
      subgroup_moved == "both" & (n_A == 1 | n_B == 1) ~ "one",
      subgroup_moved == "both" & (n_A > 1 | n_B > 1) ~ "many"
    )
  ) %>%
  group_by(event_type, split_type, group) %>%
  summarise(
    Adult_Male = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Adult_Male, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Adult_Male, 0), na.rm = TRUE),
    Subadult_Male = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Subadult_Male, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Subadult_Male, 0), na.rm = TRUE),
    Juvenile_Male = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Juvenile_Male, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Juvenile_Male, 0), na.rm = TRUE),
    Adult_Female = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Adult_Female, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Adult_Female, 0), na.rm = TRUE),
    Subadult_Female = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Subadult_Female, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Subadult_Female, 0), na.rm = TRUE),
    Juvenile_Female = sum(ifelse(subgroup_moved %in% c("A", "both"), A_Juvenile_Female, 0), na.rm = TRUE) +
      sum(ifelse(subgroup_moved %in% c("B", "both"), B_Juvenile_Female, 0), na.rm = TRUE)
  ) %>%
  pivot_longer(cols = Adult_Male:Juvenile_Female, names_to = "age_sex_class", values_to = "count")

one_many_counts$age_sex_class <- sub("_", " ", one_many_counts$age_sex_class)
one_many_counts$age_sex_class <- as.factor(one_many_counts$age_sex_class)
one_many_counts$age_sex_class <- factor(one_many_counts$age_sex_class, levels = c("Adult Female", "Adult Male", "Subadult Female", "Subadult Male", "Juvenile Female", "Juvenile Male"))

g2 <- ggplot(one_many_counts, aes(x = split_type, y = count, fill = age_sex_class)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_classic(base_size = 20) +
  theme(
    axis.text.x = ggtext::element_markdown(margin = margin(t = -20)) # Adjust margin here
  ) +
  guides(fill=guide_legend(title="Age/Sex class"))+
  scale_fill_manual(values = custom_colors) + 
  xlab(paste0(event_type, " type")) + ylab("Count") +
  scale_x_discrete(labels = function(x) glue::glue("<img src='C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/arrows/{x}.png' height='30' />")) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = ggtext::element_markdown(margin = margin(t = -10)), # Adjust margin here
    axis.text.y = element_text(size = 20), # Adjust y-axis text if needed
    axis.title.x = element_text(size = 20, margin = margin(t = 10)),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    panel.border = element_blank(),       # Remove the outline around the plot
    panel.grid.major = element_blank(),   # Remove major grid lines if desired
    panel.grid.minor = element_blank(),   # Remove minor grid lines if desired
    panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
    plot.background = element_rect(fill = "white", color = NA)   # Ensure plot background is blank
  ) +
  guides(fill = guide_legend(title = "Age/Sex class")) +
  facet_wrap(~group)

g2

ggsave(g2, file = paste0(plot_dir, group,"_", event_type, "_agesex_one_manymove_nomales.png"), width = 11, height = 8)

