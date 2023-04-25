#this script is for comparing the fission fusion events in the automated detection to the manual labels

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

group <- 'galaxy' #subdirectory where the group data is stored

codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
#groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(codedir)

#read in events
events <- read.csv(paste0('Split_mechanics/',group,'_manual_split_merge_clean.csv'), sep=';')


























#detect_fissions_and_fusions <- function(R_inner, R_outer, xs = xs, ys = ys, ts = ts, coati_ids = coati_ids, verbose = T)





