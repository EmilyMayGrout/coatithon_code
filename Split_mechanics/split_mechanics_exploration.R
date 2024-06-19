#Exploring fission-fusion mechanics
#defining the different types of splits

#LIBRARY
library(lubridate)
library(scales)
library(ggplot2)
library(patchwork)
library(dplyr)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS

#who is using (ari or emily)
user <- 'emily'

#which group - galaxy or presedente
group <- 'presedente' #subdirectory where the group data is stored

#whether to identify splits and merges automatically (if F) or use manually identified events (if T)
use_manual_events <- F

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
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
    
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
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

hist(log(fiss_speed), breaks = 100, freq = FALSE, main = " ", xlab = "Log speed during displacement (m/s)", ylab = "Density", col = "darkslategray3")


#decided to cut-off the "non-moving" from the "moving" group at 10 m based on visual inspection of the plots and from biological reasoning 
#Create column 'subgroup_move' using ifelse

dist_thresh <- 10
type <- "fission"

events$subgroup_move <- ifelse(events$A_during_disp > dist_thresh & events$B_during_disp > dist_thresh, "Both moved",
                               ifelse(events$A_during_disp > dist_thresh, "Either moved",
                                      ifelse(events$B_during_disp > dist_thresh, "Either moved",
                                             "Neither moved")))

# Create the new column to see if the distances moved are similar 
events$speed_comparison <- ifelse(abs(events$A_during_disp - events$B_during_disp) < dist_thresh, "Similar distance travelled", "Different distance travelled")

#was it a group or an individual who moved
#not sure this makes sense as it says singleton if just one of the groups contained a singleton
#events$individual_movement <- ifelse((events$A_during_disp > dist_thresh & events$n_A == 1) | 
#                                       (events$B_during_disp > dist_thresh & events$n_B == 1), 
#                                     "singleton", "multiple individuals")


if(group == 'presedente'){
  high_col <- "mediumpurple4"
} else if(group == "galaxy"){
  high_col <- "darkcyan"
}



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
      n_A == 1 & n_B == 1 ~ "1/1",
      n_A > 1 & n_B == 1 ~ "1/many",
      n_A == 1 & n_B > 1 ~ "1/many",
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
save(detailed_events, file = paste0(groupdir, group, "_detailed_events.RData"))
}else if(group == "presedente"){
  save(detailed_events, file = paste0(groupdir, group, "_detailed_events.RData"))
}

#-------------------------------------------------------
#-------------------------------------------------------

#define type of event to look at
event_type <- "fission"

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
  labs(title = paste(event_type, "heatmap"),
       x = "Subgroup Composition",
       y = "Event Type",
       fill = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
all

file_path <- file.path(plot_dir, paste(event_type, "matrix.png"))
ggsave(file_path, all, width = 10, height = 8)




#redo matrix with and without males to see how this affects the patterns observed

#find the events where males are involved 
if(group == 'presedente'){
 male_coatis <- c("Sam", "Ken", "Gen", "Lul")
} else if(group == "galaxy"){
  male_coatis <- "Gus"
}

# Function to check if all individuals in a group are males
all_males <- function(group, male_names) {
  individuals <- strsplit(group, ",\\s*")[[1]]
  all(individuals %in% male_names)
}

# Filter the dataframe
filtered_event_type <- event_type_df %>%
  rowwise() %>%
  filter(all_males(group_A, male_coatis) | all_males(group_B, male_coatis))

# Filter the dataframe for events where not all individuals leaving are males
filtered_event_type_no_males <- event_type_df %>%
  rowwise() %>%
  filter(!(all_males(group_A, male_coatis) | all_males(group_B, male_coatis)))

#choose whether want the male events or non-male events
with_males <- F

# Create a new dataframe based on the choice
df <- if (with_males) {
  filtered_event_type # DataFrame with male events
  
} else {
  filtered_event_type_no_males  # DataFrame without male events
} 


if(with_males){
  name <- "with males"
}else{
  name <- "without males"
}


contingency_table <- table(df$split_type, df$subgroup_comp)
contingency_matrix <- as.matrix(contingency_table)
# Convert the matrix to a data frame for ggplot2
contingency_df <- as.data.frame(as.table(contingency_matrix))


g <- ggplot(contingency_df, aes(Var2, Var1, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = high_col) +
  labs(title = paste(event_type, "heatmap", name),
       x = "Subgroup Composition",
       y = "Event Type",
       fill = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g


file_path <- file.path(plot_dir, paste(name, event_type, "events.png"))
ggsave(file_path, g, width = 10, height = 8)




#now want to look at the individual who moves to see if this is the singleton or the group 
#for events where both are still and one moved, was it the single individual who moved?

#choose whether to look at individual counts with or without the males 
df_filt <- filtered_event_type #filtered_fission_df or filtered_fission_df_no_males

df_filt$singleton_move <- ifelse(df_filt$split_type == "bothstill_onemove" & df_filt$subgroup_comp == "1/many" & df_filt$n_A == 1 & df_filt$A_during_disp > 10, "singleton_move",
                           ifelse(df_filt$split_type == "bothstill_onemove" & df_filt$subgroup_comp == "1/many" & df_filt$n_B == 1 & df_filt$B_during_disp > 10, "singleton_move",
                           ifelse(df_filt$split_type == "bothmove_onemove" & df_filt$subgroup_comp == "1/many" & df_filt$n_A == 1 & df_filt$A_during_disp < 10, "singleton_stop",
                           ifelse(df_filt$split_type == "bothmove_onemove" & df_filt$subgroup_comp == "1/many" & df_filt$n_B == 1 & df_filt$B_during_disp < 10, "singleton_stop",  "other"))))

df_filt <- df_filt %>% 
  filter(singleton_move != "other")

table(df_filt$singleton_move)
barplot(table(df_filt$singleton_move), col = "skyblue" )
#In Presidente group from 44 events where there is one individual involved, if the group was still, and one group moves, 38/44 events are a single individual leaving (34 of these are males). If the group were moving and one stops, 6/16 are a single individual (4 are males)

#so seems like the majority of fissions are driven by the males














### FISSION PLOTTING ###
fis <- events$event_type == "fission"

#combine the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp to make abline
combined_1 <- data.frame(events$A_subgroup_size[fis])
colnames(combined_1) <- "sub_size"
disp_1 <- data.frame(events$A_during_disp[fis])
colnames(disp_1) <- "disp"
combined_2 <- data.frame(events$B_subgroup_size[fis])
colnames(combined_2) <- "sub_size"
disp_2 <- data.frame(events$B_during_disp[fis])
colnames(disp_2) <- "disp"

combined <- rbind(combined_1,combined_2)
disp <- rbind(disp_1, disp_2)
comb <- cbind(combined, disp)


png(height = 800, width = 2500, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfission.png'))
par(mfrow=c(1,3), mar = c(10,10,10,10),(mgp=c(3,5,1))) #bottom, left, top, right)
#plot 1:
plot(events$A_during_disp[fis], events$B_during_disp[fis], pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Distance traveled during fission event (m)', cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "aquamarine3",mgp=c(5,2,.5))
abline(lm(events$B_during_disp[fis] ~ events$A_during_disp[fis]))

#plot 2:
#combined the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp so I could add the abline
plot(comb$disp, comb$sub_size, pch = 20, ylim = c(0, 10), xlab = "Distance traveled during fission event (m)", ylab = "Sub-group size", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2",mgp=c(5,2,.5))
#abline(lm(comb$sub_size ~ comb$disp))

hist(events$AB_before_disp[fis], main = "Distance traveled 10 minutes before displacement (m)", xlab = '', col = "aquamarine3",cex.axis = 3, cex.lab = 3, cex.main = 3,mgp=c(5,2,.5))

#filtering split angle times when distance of full group travelled more that 20m
#events_dist_20 <- events[events$AB_before_disp > 20,]
#hist(events_dist_20$split_angle[events_dist_20$event_type == "fission"], main = "Split angle between sub-groups after fission (degrees)", xlab = '', cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral3",mgp=c(5,2,.5))

dev.off()



#filter dataframe to columns needed for the distance travelled and order the subgroup by size (so the larger subgroup is on the X axis and the smaller subgroup is on the Y axis)

fission_df <- data.frame(cbind(events$A_during_disp[fis], events$B_during_disp[fis], events$A_subgroup_size[fis], events$B_subgroup_size[fis], events$AB_before_disp[fis]))
colnames(fission_df) <- c("A_during_disp", "B_during_disp", "A_subgroup_size", "B_subgroup_size", "AB_before_disp")

#create column for the distnace travelled for the larger subgroup
fission_df$Distance_Larger_Group <- ifelse(fission_df$A_subgroup_size >= fission_df$B_subgroup_size, fission_df$A_during_disp, fission_df$B_during_disp)
#create column for the distance travelled of the smaller subgroup
fission_df$Distance_Smaller_Group <- ifelse(fission_df$A_subgroup_size >= fission_df$B_subgroup_size,fission_df$B_during_disp, fission_df$A_during_disp)

png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfission_groupsizeorder_lm.png'))
par(mar = c(8,8,2,2))
plot(fission_df$Distance_Larger_Group, fission_df$Distance_Smaller_Group, pch = 20, xlab = "Distance travelled by larger subgroup (m)", ylab = "Distance travelled by smaller subgroup (m)", main = '', cex = 4, cex.axis = 2.5, cex.lab = 2, col = "aquamarine3",mgp=c(5,2,0))
abline(lm(fission_df$Distance_Smaller_Group ~ fission_df$Distance_Larger_Group))
dev.off()

png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_beforefission_groupsizeorder.png'))
par(mar = c(8,8,2,2), mgp = c(3,0,-1.7))
hist(fission_df$AB_before_disp, main = '', xlab = "Distance traveled 10 minutes before fission (m)", ylab = '', col = "aquamarine4",cex.axis = 2.5, cex.lab = 2, cex.main = 3, breaks = 25)
dev.off()



### PLOTTING EACH AGE CLASS IN ONE GRAPH ###
#changing the xlim for Galaxy to 60 (Presedente is 80)
a <- 0.5
png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_ages_duringfission.png'))
par(mfrow=c(1,1), mar = c(10,10,2,2),(mgp=c(3,5,1))) #bottom, left, top, right)
#combine the 3 plots above into one plot
plot(events$A_during_disp[fis], events$A_proportion_adults[fis], pch = 19, ylab = "Proportion of age class in sub-group", xlab="Distance traveled during fission event (m)",xlim = c(0, 80), ylim = c(0,1), cex = 4, cex.axis = 2, cex.lab = 2, cex.main = 3, col = alpha("hotpink4",a), mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_proportion_adults[fis], pch = 19,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("hotpink4",a))
#plot distance traveled depending on proportion of sub-adults in subgroup
points(events$A_during_disp[fis], events$A_proportion_subadults[fis], pch = 17,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("orange3", a),mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_proportion_subadults[fis], pch = 17,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("orange3", a))

#plot distance traveled depending on proportion of juveniles in subgroup
points(events$A_during_disp[fis], events$A_proportion_juveniles[fis], pch = 15, ylim = c(0,1), cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("cadetblue3",a), mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_proportion_juveniles[fis], pch = 15,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("cadetblue3", a))
legend(x = "topright",         
       legend = c("Adult", "Sub-adult", "Juvenile"),
       pch = c(19, 17, 15),
       col = c("hotpink4", "orange3", "cadetblue3"),
       cex = 2, pt.cex = 2,
       bty = "n")
dev.off()



#removed plot as the turn_angle gives less info than the split angle (because not positive or negative)
plot(events$turn_angle_A[fis], events$turn_angle_B[fis], pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Turn angle after fission event (degrees)',col = "coral3",mgp=c(5,2,.5))

#same as plot 2 from above for loop, but the data here is directly from the events dataframe (rather than my rbind data frame to make the abline)
plot(events$A_during_disp[fis],events$A_subgroup_size[fis], pch = 20, ylim = c(0, 10), xlab = "Distance traveled during fission event (m)", ylab = "Sub-group size", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2",mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_subgroup_size[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2")
abline(lm(events$B_subgroup_size[fis] ~ events$B_during_disp[fis]))


#plotting sub-group age and distance traveled during event
#combining the A group and the B group as it's not about the group id, but the average age (or number of adults) in those groups
plot(events$A_during_disp[fis], events$A_average_grp_age[fis], xlim = c(0,51), ylim = c(2,3), pch = 20, ylab = "Group average age class (1 = juvenile, 2 = subadult, 3 = adult)", xlab = "Distance traveled during fission event (m)")
points(events$B_during_disp[fis], events$B_average_grp_age[fis], pch = 20)



png(height = 600, width = 2400, units = 'px', filename = paste0(plot_dir,'inds_dist_traveled_duringfission.png'))

par(mfrow=c(1,4), mar = c(8,8,5,8),(mgp=c(3,3,1)))
#plot distance traveled depending on number of adults in subgroup
plot(events$A_during_disp[fis], events$A_n_adults[fis], pch = 20, xlab = "Distance traveled during fission event (m)", ylab = "Number of adults in subgroup", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_n_adults[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3)

#plot distance traveled depending on number of subadults in subgroup
plot(events$A_during_disp[fis], events$A_n_subadults[fis], pch = 20,  xlab = "Distance traveled during fission event (m)", ylab = "Number of subadults in subgroup", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_n_subadults[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3)

#plot distance traveled depending on number of juveniles in subgroup
plot(events$A_during_disp[fis], events$A_n_juveniles[fis], pch = 20, xlab = "Distance traveled during fission event (m)", ylab = "Number of juveniles in subgroup", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_n_juveniles[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3)

#change ylim to 10 when Galaxy group, 20 when Presedente
#plot distance traveled depending on number of individuals in subgroup
plot(events$A_during_disp[fis],events$A_subgroup_size[fis], pch = 20, xlab = "Distance traveled during fission event (m)", ylab = "Number of individuals in subgroup", ylim = c(0, 20), cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, mgp=c(5,2,.5), col = "cyan3")
points(events$B_during_disp[fis], events$B_subgroup_size[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = "cyan3")

dev.off()

#plot distance traveled depending on proportion of adults in subgroup
plot(events$A_during_disp[fis], events$A_proportion_adults[fis], pch = 20)
points(events$B_during_disp[fis], events$B_proportion_adults[fis], pch = 20)


# ----------------------------------------------------------------------
### FUSION PLOTTING ###

fus <- events$event_type == "fusion"

#need to think about what sort of plots to do...

#combine the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp to make abline
combined_fus_1 <- data.frame(events$A_subgroup_size[fus])
colnames(combined_fus_1) <- "sub_size"
disp_fus_1 <- data.frame(events$A_during_disp[fus])
colnames(disp_fus_1) <- "disp"
combined_fus_2 <- data.frame(events$B_subgroup_size[fus])
colnames(combined_fus_2) <- "sub_size"
disp_fus_2 <- data.frame(events$B_during_disp[fus])
colnames(disp_fus_2) <- "disp"

combined_fus <- rbind(combined_fus_1,combined_fus_2)
disp_fus <- rbind(disp_fus_1, disp_fus_2)
comb_fus <- cbind(combined_fus, disp_fus)

##want to look at distance traveled before a fusion 
png(height = 1000, width = 2000, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfusion.png'))
par(mfrow=c(1,2), mar = c(10,10,10,10),(mgp=c(3,5,1))) #bottom, left, top, right)

#plot 1:
plot(events$A_during_disp[fus], events$B_during_disp[fus], pch = 20, xlab = "", ylab = "", main = 'Distance traveled during fusion event (m)', col = "cadetblue4",cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,mgp=c(5,2,.5))
abline(lm(events$B_during_disp[fus] ~ events$A_during_disp[fus]))

#plot 2:
#combined the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp so I could add the abline
plot(comb_fus$disp, comb_fus$sub_size, pch = 20, ylim = c(0, 10), xlab = "Distance traveled during fusion event (m)", ylab = "Sub-group size", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "hotpink4",mgp=c(5,2,.5))
#abline(lm(comb_fus$sub_size ~ comb_fus$disp))

#hist(events$split_angle[fus], xlab = "Angle between sub-groups during fusion (degrees)", main = " ",cex.axis = 3, cex.lab = 3, cex.main = 3,col = "hotpink4",mgp=c(5,2,.5))

dev.off()


#filter dataframe to columns needed for the distance travelled and order the subgroup by size (so the larger subgroup is on the X axis and the smaller subgroup is on the Y axis)

fusion_df <- data.frame(cbind(events$A_during_disp[fus], events$B_during_disp[fus], events$A_subgroup_size[fus], events$B_subgroup_size[fus]))
colnames(fusion_df) <- c("A_during_disp", "B_during_disp", "A_subgroup_size", "B_subgroup_size")

#create column for the distance travelled for the larger subgroup
fusion_df$Distance_Larger_Group <- ifelse(fusion_df$A_subgroup_size >= fusion_df$B_subgroup_size, fusion_df$A_during_disp, fusion_df$B_during_disp)
#create column for the distance travelled of the smaller subgroup
fusion_df$Distance_Smaller_Group <- ifelse(fusion_df$A_subgroup_size >= fusion_df$B_subgroup_size,fusion_df$B_during_disp, fusion_df$A_during_disp)


png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfusion_groupsizeorder.png'))
par(mar = c(8,8,2,2))
plot(fusion_df$Distance_Larger_Group, fusion_df$Distance_Smaller_Group, pch = 20, xlab = "Distance travelled by larger subgroup (m)", ylab = "Distance travelled by smaller subgroup (m)", main = '', cex = 4, cex.axis = 2.5, cex.lab = 2, col = "cadetblue3",mgp=c(5,2,0))
#abline(lm(fusion_df$Distance_Smaller_Group ~ fusion_df$Distance_Larger_Group))
dev.off()

#for changing the alpha values
library("scales") 

#change xlim for galaxy to 60, presedente is 100

a <- 0.5
png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_ages_duringfusion.png'))
par(mfrow=c(1,1), mar = c(10,10,2,2),(mgp=c(3,5,1))) #bottom, left, top, right)
#combine the 3 plots above into one plot
plot(events$A_during_disp[fus], events$A_proportion_adults[fus], pch = 19, ylab = "Proportion of age class in sub-group", xlab="Distance traveled during fusion event (m)",xlim = c(0, 100), ylim = c(0,1), cex = 4, cex.axis = 2, cex.lab = 2, cex.main = 3, col = alpha("hotpink4",a), mgp=c(5,2,.5))
points(events$B_during_disp[fus], events$B_proportion_adults[fus], pch = 19,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("hotpink4",a))
#plot distance traveled depending on proportion of sub-adults in subgroup
points(events$A_during_disp[fus], events$A_proportion_subadults[fus], pch = 17,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("orange3", a),mgp=c(5,2,.5))
points(events$B_during_disp[fus], events$B_proportion_subadults[fus], pch = 17,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("orange3", a))

#plot distance traveled depending on proportion of juveniles in subgroup
points(events$A_during_disp[fus], events$A_proportion_juveniles[fus], pch = 15, ylim = c(0,1), cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("cadetblue3",a), mgp=c(5,2,.5))
points(events$B_during_disp[fus], events$B_proportion_juveniles[fus], pch = 15,cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3, col = alpha("cadetblue3", a))
legend(x = "topright",         
       legend = c("Adult", "Sub-adult", "Juvenile"),
       pch = c(19, 17, 15),
       col = c("hotpink4", "orange3", "cadetblue3"),
       cex = 2, pt.cex = 2,
       bty = "n")
dev.off()

#the plot above isn't a good way of looking at the age class and distance traveled 
#now want to plot distance traveled with the proportion of age calls based on the total number of that age class
events$total_n_adults <- events$A_n_adults + events$B_n_adults
events$total_n_subadults <- events$A_n_subadults + events$B_n_subadults
events$total_n_juveniles <- events$A_n_juveniles + events$B_n_juveniles

#filter events to remove instances when one individual is fissioning or fusioning
events_groups <- subset(events, events$A_subgroup_size > 1 & events$B_subgroup_size > 1, drop = TRUE)

#filter events_groups to when one group moves and the other stays (moved less than 15m during event)
events_groups <- subset(events_groups, events_groups$A_during_disp < 20 | events_groups$B_during_disp  < 20, drop = TRUE)


#adults
plot(events_groups$A_during_disp[fis], (events_groups$A_n_adults[fis]/events_groups$total_n_adults[fis]), xlim = c(0,90))
points(events_groups$B_during_disp[fis], (events_groups$B_n_adults[fis]/events_groups$total_n_adults[fis]))
#subadults
plot(events_groups$A_during_disp[fis], (events_groups$A_n_subadults[fis]/events_groups$total_n_subadults[fis]), xlim = c(0,90))
points(events_groups$B_during_disp[fis], (events_groups$B_n_subadults[fis]/events_groups$total_n_subadults[fis]))

