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
group <- 'galaxy' #subdirectory where the group data is stored

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
  ylims = c(0, 0.06)
  xlims = c(0, 80)
} else if(group == "galaxy"){
  high_col <- "darkolivegreen4"
  comp_col <- "darkolivegreen1"
  break_1 <- 10
  break_2 <- 40
  break_3 <- 5
  break_4 <- 20
  break_5 <- 50
  break_6 <- 50
  ylims = c(0, 0.11)
  xlims = c(0, 60)
  
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


events_clean <- events[events$n_A != 1 & events$n_B != 1, ]
events_clean$length_during <- events_clean$end_time - events_clean$start_time

mean(na.omit(events_clean$length_during))
sd(na.omit(events_clean$length_during))


events_clean <- events_clean[is.finite(events_clean$length_during), ]

ggplot(events_clean, aes(x = length_during, fill = event_type)) +
  geom_histogram(bins = break_1, color = "white", position = "dodge") +
  labs(
    x = "During duration (seconds)",
    y = "Frequency",
    fill = "Event type" 
  ) +
  theme_classic() +                 
  theme( 
    axis.text = element_text(size = 20),  
    axis.title = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16) 
  ) +
  scale_fill_manual(values = c("fission" = high_col, "fusion" = comp_col))



ggsave(paste0(plot_dir, "during_duration.png"), width = 8, height = 5)

