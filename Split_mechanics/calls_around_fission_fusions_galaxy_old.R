#this script is to look at the calls used and call rates before and after events

#read in labels 

wd <- "C:/Users/egrout/Dropbox/coaticalls/Galaxy_labels/completed_labels/labels_cleaned_25.02.24/"
plot_dir <- "C:/Users/egrout/Dropbox/coaticalls/results/"

setwd <- wd

#read in events - RData was made in split_mechanics_exploration 
load('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxy_manual_events_withinfo2.RData')

#create column for the distnace travelled for the larger subgroup
events$Distance_Larger_Group <- ifelse(events$A_subgroup_size >= events$B_subgroup_size, events$A_during_disp, events$B_during_disp)
#create column for the distance travelled of the smaller subgroup
events$Distance_Smaller_Group <- ifelse(events$A_subgroup_size >= events$B_subgroup_size,events$B_during_disp, events$A_during_disp)

#read in coati IDs 
load('C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_coati_ids.RData')


#LIBRARIES
library(stringr)
library(hms)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(gridExtra)

# get a list of all the CSV files in the folder
files <- list.files(wd, pattern = "*.csv")

#create data frame with 0 rows and 3 columns
all_data <- data.frame(matrix(ncol = 8, nrow = 0))
#provide column names
colnames(all_data) <- c("label","Start","Duration","Time","Format","Type","Description","file_name")

# loop through each CSV file and add it to the dataframe
#i = 67

for (i in 1:length(files)) {
  # read in the CSV data as a tibble
  # using header = TRUE assumes the first row of each CSV file is a header with column names
  file_data <- read.csv(paste0(wd, files[i]), header = T)
  
  # add a column with the row names (i.e. the name of the CSV file)
  file_data$file_name <- files[i]
  #only keeping necessary info of the file name
  file_data$file <- str_sub(file_data$file_name,end=11)
  file_data$date <- str_sub(file_data$file_name, start = 13, end = 20)
  colnames(file_data)[colnames(file_data) == "Name"] <- "label"
  
  # add the data to the all_data dataframe
  all_data <- rbind(all_data, file_data)
}


#remove rows which contain the date 
all_data <- all_data[!grepl(":", all_data$label),]

#remove file_name column 
all_data <- all_data[,-7]

#add ID column
all_data$id <- str_sub(all_data$file, end=5)

#make date colum in POSIXct
all_data$date <- as.Date(sub("(..)$", "20\\1", all_data$date), "%d.%m.%Y")

#make time column
all_data$Start <- all_data$Start
table(str_length(all_data$Start))

#because the length of the time column is different due to some times less than an hour, need to split the data to get the times and then rbind them
t <- filter(all_data, nchar(Start) == 8)
t$time <- paste0("6:", t$Start)
t$time <- as_hms(t$time)

s <- filter(all_data, nchar(Start) == 9)
s$time <- paste0("6:", s$Start)
s$time <- as_hms(s$time)

v <- filter(all_data, nchar(Start) == 11)
v$time <- as_hms(as_hms(v$Start) + as_hms('6:00:00.000'))

all_data_hms <- rbind(t,s,v)

#add as.POSIXct
all_data_hms$datetime <- as.POSIXct(paste(all_data_hms$date, all_data_hms$time), format = "%Y-%m-%d %H:%M:%OS")

#write.csv(all_data_hms, "C:/Users/egrout/Dropbox/coaticalls/processed/all_data_hms.csv", row.names = F)


#in case want to split by id for plotting
#f1 <- all_data_hms %>% group_split(id)

#for now, to solve the time synch (though this needs correcting properly later) I'm adding 1 minute to all calls

all_data_hms$datetime_new <- ifelse(all_data_hms$label %in% c('fission', 'fusion'), all_data_hms$datetime,all_data_hms$datetime + lubridate::seconds(60))

all_data_hms$datetime_new <- as.POSIXct(all_data_hms$datetime_new)

#combining the contact calls 
all_data_hms$label[all_data_hms$label == "chirpgr"] <- "contact call"
all_data_hms$label[all_data_hms$label == "chirp grunt"] <- "contact call"
all_data_hms$label[all_data_hms$label == "chirp click"] <- "contact call"
all_data_hms$label[all_data_hms$label == "click grunt"] <- "contact call"
all_data_hms$label[all_data_hms$label == "click"] <- "contact call"
all_data_hms$label[all_data_hms$label == "chirp"] <- "contact call"


#combine aggressive calls
all_data_hms$label[all_data_hms$label == "chitter"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "squeal"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "squeal chitter"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "squeal chitter x"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "squeal chitters"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "low squeal"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "chitter x"] <- "aggression call"
all_data_hms$label[all_data_hms$label == "squeal chittering"] <- "aggression call"

#add coati names to column based on IDs
all_data_hms$name[all_data_hms$id == "G9463"] <- "Estrella"
all_data_hms$name[all_data_hms$id == "G9476"] <- "Luna"
all_data_hms$name[all_data_hms$id == "G9474"] <- "Saturno"
all_data_hms$name[all_data_hms$id == "G9464"] <- "Venus"
all_data_hms$name[all_data_hms$id == "G9470"] <- "Orbita"
all_data_hms$name[all_data_hms$id == "G9475"] <- "Pluto"
all_data_hms$name[all_data_hms$id == "G9480"] <- "Cometa"
all_data_hms$name[all_data_hms$id == "G9466"] <- "Lucero"
all_data_hms$name[all_data_hms$id == "G9471"] <- "Planeta"
all_data_hms$name[all_data_hms$id == "G9460"] <- "Quasar"
all_data_hms$name[all_data_hms$id == "G9467"] <- "Gus"


#change label time for correct time!! (I made a mistake when asking the HiWis to label around the fusion event, I said at 07:27 when it should have been at 07:57)
all_data_hms$datetime_new[all_data_hms$datetime == as.POSIXct("2021-12-27 07:27:00")] <- as.POSIXct("2021-12-27 07:57:00")

#BUT because the results of this event are incorrect as not all the time around the event have been labelled, this fusion event needs to be removed from the analysis

all_data_hms <- all_data_hms[!(all_data_hms$datetime_new == as.POSIXct("2021-12-27 07:57:00")),]


#create a data frame to hold output data: columns are fission time, ind id, chirp grunt rate before, chirpgr rate after, chitter rate before, chitter after, etc.
#df_out <- data.frame(time = NA, ind = NA, fis_or_fus = NA, chirp_bef = NA, chirp_aft = NA, chirpgr_bef = NA, chirpgr_aft = NA, chitter_bef = NA, chitter_aft = NA) #fake first row might be necessary

df_out <- data.frame()
# i = 1
# j = 2
# n = 5

for (i in unique(all_data_hms$id)){

  #extract each individuals data for all days
  ind_i <- all_data_hms[all_data_hms$id == i,]
  
  #for loop through each day
  for (j in unique(ind_i$date)){
    day_j <- ind_i[ind_i$date == as.Date(j),]
    
    for (n in 1:nrow(day_j)){
      #go through each row of dataframe, if the label is "fission", then extract all labels 15 mins before and after that row number to get a table of call counts
      row_n <- day_j[n,]
      label <- row_n$label
      if(label %in% c("fission", 'fusion')){ #if row_n is either a fission or a fusion
        
        #if the row is a fission, get the time of that fission, then get the time of calls 15 mins before and 15 mins after
        fis_time <- day_j$datetime_new[n]
        time <- lubridate::seconds(900) #15 mins
        first_time <- fis_time - time
        last_time <-  fis_time + time
        
        labels_bef <- day_j[day_j$datetime_new >= first_time & day_j$datetime_new <= fis_time,]
        labels_aft <-day_j[day_j$datetime_new >= fis_time & day_j$datetime_new <= last_time,]
        
        #get the number of each type of call before
        label_count_bef <- data.frame(table(labels_bef$label))
        label_count_bef$rate <- (label_count_bef$Freq)/(as.numeric(time)/60) #rate in mins
      
        #get the number of each type of call after
        label_count_aft <- data.frame(table(labels_aft$label))
        label_count_aft$rate <- (label_count_aft$Freq)/(as.numeric(time)/60) #rate in mins
        
        #create a row of your data frame above and append it using rbind
        row <- data.frame(time = fis_time, ind = i, fis_or_fus = row_n, contact_call_aft = NA, contact_call_bef = NA, aggression_call_aft = NA, aggression_call_bef = NA)
        
       #for loop to extract the call rates for the calls we're interested in:
        for(f_call in c("contact call", "aggression call")){
          rate_aft <- label_count_aft$rate[label_count_aft$Var1==f_call]
          rate_bef <- label_count_bef$rate[label_count_bef$Var1==f_call]
        
          row[, paste0(f_call, "_aft")] <- rate_aft
          row[, paste0(f_call, "_bef")] <- rate_bef
      }
        df_out <- rbind(df_out, row) #bind the row to the current data frame
        
      }
        else{NULL} #not needed probably
  }
 }
}
    


#for plotting will pivot this table longer and remove some of the columns
df_out_filt <- df_out[,c(1,2,3,15,20:23)]

#pivot longer
df_out_long <- df_out_filt %>%
  pivot_longer(!c(time, ind, fis_or_fus.label, fis_or_fus.name), names_to = "call", values_to = "rate")
#replace NA with 0
df_out_long[is.na(df_out_long)] <- 0

#make column for before and after
df_out_long[c('call', 'bef_aft')] <- str_split_fixed(df_out_long$call, '_', 2)

#write name of bef and aft longer
df_out_long$bef_aft[df_out_long$bef_aft == "bef"] <- "Before"
df_out_long$bef_aft[df_out_long$bef_aft == "aft"] <- "After"


#change factor levels so before is shown before after in plot
df_out_long$bef_aft <- factor(df_out_long$bef_aft, levels = c("Before", "After"))

g <- ggplot(data = df_out_long[df_out_long$fis_or_fus.label == "fusion",], aes(x = call, y = rate))+
  geom_boxplot(aes(fill = bef_aft))+ 
  geom_point(position=position_jitterdodge(), size = 0.5, color = "gray3", aes(fill=bef_aft, alpha = 0.8))+
  scale_fill_manual(values=c("indianred1", "indianred4"))+
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_text(size=14), axis.title = element_text(size = 14), legend.text = element_text(size = 14))+ 
  xlab(" ") +
  ylab("Call rate (per minute)")+
  scale_alpha(guide = 'none')#+
  facet_wrap(~time)

g


#for fissions: "darkolivegreen2", "darkolivegreen"
#for fusions: "indianred1", "indianred4"

ggsave(paste0(plot_dir, "fusion_perevent_15mins.png"), width = 50, height = 20, units = "cm", g)





#look at one event to see who is vocal, then see if they moved or not
one_event <- df_out_long[df_out_long$fis_or_fus.label == "fission",]

time <- "2021-12-28 06:55:00 CET"
#time <- "2022-01-04 07:02:00 CET"
#time <- "2022-01-02 06:43:00 CET"
#time <- "2021-12-28 07:29:00 CET"

g <- ggplot(data = one_event[one_event$time == time,], aes(x = call, y = rate, color = bef_aft))+ 
  geom_point(size = 5) +
  scale_colour_manual(breaks = one_event[one_event$time == time,]$bef_aft, values = c("darkolivegreen2", "darkolivegreen"))+
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text=element_text(size=14), axis.title = element_text(size = 14), legend.text = element_text(size = 14))+ 
  xlab(" ") +
  ylab("Call rate (per minute)")+
  scale_alpha(guide = 'none')+
  facet_wrap(~fis_or_fus.name)

g

#for fissions: "darkolivegreen2", "darkolivegreen"
#for fusions: "indianred1", "indianred4"


ggsave(paste0(plot_dir, "fission_0655_perind_15mins.png"), width = 35, height = 20, units = "cm", g)





#---------------------------------------------------------------------------------------------------

#look at call rates of individuals who are the "leavers"

#need to put events df and df_out_filt df on same time zone
df_out_filt$UTC_time <- df_out_filt$time + 5*60*60 # add 5 hours
#filter for events that we have call rates results for
events_filt <- events[events$datetime %in% unique(df_out_filt$UTC_time), ]

#colnames(events_filt)
#remove some column that are unnecessary 
events_filt <- events_filt[,c(2,3,4,5,6,9,10,12,14,15,20,21,22)]

#split calling behaviour by event
#f1 <- df_out_filt %>% group_split(time)

#look at one event for the calling rate of each group member and see if there's a relationship between which subgroup they were in
#look at the first event: 28.12.21 11:55:00


events_filt_all <- data.frame()

for (i in 1:nrow(events_filt)){

  #get the calling behaviour of the first event
  call_df_event_1 <- df_out_filt[df_out_filt$UTC_time == unique(df_out_filt$UTC_time)[i],]

  #get the events movement info by finding the same date of the calling event to the events_filt dataframe
  events_filt_1 <- events_filt[events_filt$datetime %in% unique(df_out_filt$UTC_time)[[i]],]

  #I want to put each individual in group A and B into a dataframe where each ind is a row and they have another column for distance travelled and subgroup ID

  #extract the IDs from the list for event i to each row of a dataframe
  id_A <- data.frame(matrix(unlist(events_filt_1$group_A_idxs), ncol = length(events_filt_1$group_A_idxs)))
  #add column for distance travelled ect.
  id_A <- cbind(id_A, events_filt_1$n_A, events_filt_1$A_during_disp, events_filt_1$AB_before_disp, events_filt_1$datetime)
  colnames(id_A) <- c("id", "subgroup_size", "during_dist", "before_dist", "datetime")
  #do the same for group B
  id_B <- data.frame(matrix(unlist(events_filt_1$group_B_idxs), ncol = length(events_filt_1$group_B_idxs)))
  colnames(id_B) <- "id"
  id_B <- cbind(id_B, events_filt_1$n_B, events_filt_1$B_during_disp, events_filt_1$AB_before_disp, events_filt_1$datetime)
  colnames(id_B) <- c("id", "subgroup_size", "during_dist", "before_dist", "datetime")

  #rbind the dataframes

  event_1 <- rbind(id_A, id_B)
  event_1$id_indx <- event_1$id

  #get the code for the individuals in the id column
  coati_ids$id_indx <- 1:nrow(coati_ids)

  #merge the ID column with the event_1 column 
  event_1_merge <- merge(event_1, coati_ids, by="id_indx")
  #add G to tag ID
  event_1_merge$ind <- paste0("G", event_1_merge$tag_id)

  #merge to call_df_event_1
  event_calls_merge <- merge(call_df_event_1, event_1_merge, by = "ind")

  #filter to columns we want
  #event_calls_merge <- event_calls_merge[,-c(8,12)]

  #put data in long format for plotting
  event_calls_merge_long <- event_calls_merge %>%
    pivot_longer(!c(time, UTC_time, datetime, ind, fis_or_fus.name,id_indx, tag_id, id, fis_or_fus.label, subgroup_size, during_dist, before_dist, name, age,sex,color), 
                 names_to = "call", values_to = "rate")
  #replace NA with 0
  event_calls_merge_long[is.na(event_calls_merge_long)] <- 0

  #make column for before and after
  event_calls_merge_long[c('call', 'bef_aft')] <- str_split_fixed(event_calls_merge_long$call, '_', 2)

  #write name of bef and aft longer
  event_calls_merge_long$bef_aft[event_calls_merge_long$bef_aft == "bef"] <- "Before"
  event_calls_merge_long$bef_aft[event_calls_merge_long$bef_aft == "aft"] <- "After"

  event_calls_merge_long$event_id <- i 
  
  #make column for which group travelled further
  thresh1 <- mean(event_calls_merge_long$during_dist)
  
  # Create a new column based on the threshold -- this probably needs to be changed to whether the short_travel group doesn't move by a certain amount (as both groups could move similar amounts)
  event_calls_merge_long$status <- ifelse(event_calls_merge_long$during_dist > thresh1, "long_travel", "short_travel")
  
  #make column for bigger and smaller subgroup
  thresh2 <- mean(event_calls_merge_long$subgroup_size)
  
  # Create a new column based on the threshold
  event_calls_merge_long$group_size <- ifelse(event_calls_merge_long$subgroup_size > thresh2, "bigger", "smaller")
  
  events_filt_all <- rbind(events_filt_all, event_calls_merge_long)

}

#change factor levels so before is shown before after in plot
events_filt_all$bef_aft <- factor(events_filt_all$bef_aft, levels = c("Before", "After"))

one_event <- events_filt_all[events_filt_all$fis_or_fus.label == "fusion",]
#one_event <- one_event[one_event$time == "2021-12-27 07:27:00 CET",]

#boxplot the call rates of the larger group by the smaller group 
g <- ggplot(data = one_event[one_event$time == "2021-12-28 07:29:00",], aes(x = call, y = rate))+
  geom_boxplot(aes(fill = bef_aft))+ 
  geom_point(position=position_jitterdodge(), size = 0.5, color = "gray3", aes(fill=bef_aft, alpha = 0.8))+
  scale_fill_manual(values=c("indianred1", "indianred4"))+
  facet_wrap(~status + group_size, ncol = 1)+
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.text=element_text(size=14), axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14))+ 
  xlab(" ") +
  ylab("Call rate (per minute)")+
  scale_alpha(guide = 'none')
g

#for fissions: "darkolivegreen2", "darkolivegreen"
#for fusions: "indianred1", "indianred4"


ggsave(paste0(plot_dir, "fusion_0729_grpsize.png"),width = 15, height = 25, units = "cm", g)

#leaving group increases call rate before moving, especially if the leaving group is bigger
#should look at which individuals are calling at the higher rate and who follows them

#made animations for these events in gganimate_galaxy.R


#looking at other interaction plots
ggplot(data = events_filt_all[events_filt_all$fis_or_fus.label == "fission",], aes(x = call, y = rate))+
  geom_boxplot(aes(fill = bef_aft))+ 
  geom_point(position=position_jitterdodge(), size = 0.5, color = "gray3", aes(fill=bef_aft, alpha = 0.8))+
  facet_wrap(~ind+group_size)

ggplot(data = events_filt_all[events_filt_all$fis_or_fus.label == "fission",], aes(x = call, y = rate))+
  geom_boxplot(aes(fill = bef_aft))+ 
  geom_point(position=position_jitterdodge(), size = 0.5, color = "gray3", aes(fill=bef_aft, alpha = 0.8))+
  facet_wrap(~ind+status)


#to see which individuals more often leave
table <- as.data.frame(table(events_filt_all$ind, events_filt_all$status))
#because each ind has 4 rows per events, to get the correct count, need to divide by 4
table$Freq <- (table$Freq)/4
colnames(table) <- c("id", "travel", "count")
#make a stacked barplot of these results
g <- ggplot(table, aes(fill=travel, x=id, y=count)) + 
  geom_bar(position="stack", stat="identity")
g


ggsave(paste0(plot_dir, "inds_travel.png"), g)








#-----------------------------------------------------------------------------------------------------

#zero augmented gamma distribution regression model 
library(brms)


hist(df_out_long$rate)

#m1 <- brm(rate ~ bef_aft + fis_or_fus.label + (1|call) + (1|ind),
          # data = df_out_long,
          # family = hurdle_gamma,
          # control = list(adapt_delta = 0.99999,
          #                max_treedepth = 14),
          # chains = 4, iter = 1000, warmup = 500,
          # cores = 4, seed = 520, 
          # backend = "cmdstanr")




#-------------------------------------------------------------------------------------------------------------------------
#look at one fission event to see call rates of all individuals involved

#fission event at 06:55 on the 28.12.21

fission_1_281221 <- all_data_hms[all_data_hms$date == "2021-12-28",]

#filter to 15 mins before and after fission event
fission_1_281221$event <- "out"
fission_1_281221$event[fission_1_281221$time > as_hms('6:40:00.000') & fission_1_281221$time < as_hms('6:55:00.000')] <- "1_bef"
fission_1_281221$event[fission_1_281221$time > as_hms('6:55:00.000') & fission_1_281221$time < as_hms('7:10:00.000')] <- "2_aft"

#filter to around event
f1 <- fission_1_281221[fission_1_281221$event == c("1_bef", "2_aft"),]
#filter to the main call types
f1 <- f1[f1$label == c("chirp", "chirpgr", "chitter"),]

ggplot(data = f1, aes(x = label, fill = event, group = event))+
  geom_bar(position = "dodge")+
  facet_wrap(~id)







