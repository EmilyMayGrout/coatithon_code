
#code to figure out average time to fix for the 1/10 min data and 1/1h data
#spiders might be 60s, coatis might be 90s
library(lubridate)
library(tidyverse)


presedente_all <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2023/presedente/Coati Presedente CCAS BCI 2023.csv")

#filter to columns we need
pres_all_cut <- presedente_all[, c(19, 36, 42)]

#get time difference between rows based on each individual seperately
pres_all_cut <- pres_all_cut %>%
  arrange(individual.local.identifier, study.local.timestamp) %>%
  group_by(individual.local.identifier) %>%
  mutate(diff = strptime(study.local.timestamp, "%Y-%m-%d %H:%M:%S") - lag(strptime(study.local.timestamp, "%Y-%m-%d %H:%M:%S"), default = strptime(study.local.timestamp, "%Y-%m-%d %H:%M:%S")[1]),diff_secs = as.numeric(diff, units = 'secs'))

#filter to the times when the diff time is greater than 1 (excluding the 1Hz period)
pres_all_cut <- pres_all_cut[pres_all_cut$diff_secs> 60 & pres_all_cut$diff_secs < 6000,]
pres_all_cut <- pres_all_cut[, c(1:3, 6)]

#distribution of the 1/10min and 1/hour rates
hist(pres_all_cut$diff_secs, breaks = 100)

#distribution of the ttf at the 1/10 mins
hist(pres_all_cut$eobs.used.time.to.get.fix[pres_all_cut$diff_secs < 1000], breaks = 100, main = "1/10 min GPS", xlab = "Time to Fix")
mean(pres_all_cut$eobs.used.time.to.get.fix[pres_all_cut$diff_secs < 1000])
#53s
min10 <- pres_all_cut$eobs.used.time.to.get.fix[pres_all_cut$diff_secs < 1000]
#What proportion of the daytime (10minutes) fixes would we have lost if the TTF had been set to 60, 90 and 120 seconds, respectively
min10_60s <- min10[min10 < 60]
length(min10_60s)/length(min10) #69%

min10_90s <- min10[min10 < 90]
length(min10_90s)/length(min10) #80%

min10_120s <- min10[min10 < 120]
length(min10_120s)/length(min10) #86%


#distribution of the ttf at the 1/hour 
hist(pres_all_cut$eobs.used.time.to.get.fix[pres_all_cut$diff_secs > 1000], breaks = 100, main = "1/hour GPS", xlab = "Time to Fix")
mean(pres_all_cut$eobs.used.time.to.get.fix[pres_all_cut$diff_secs > 1000])
#35s


