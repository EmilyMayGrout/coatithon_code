#this code is making scatter plot for subgroup sizes for both groups for each hour of the day

library(ggplot2)
library(hms)
library(dplyr)

plot_dir <- "C:/Users/egrout/Dropbox/coatithon/results/"

#data saved from fission_fusion_galaxy_V1 and fission_fusion_presedente_V1
load("C:/Users/egrout/Dropbox/coatithon/processed/n_subs_gal.Rdata")
load("C:/Users/egrout/Dropbox/coatithon/processed/n_subs_pres.Rdata")

n_subs_gal$group <- "Galaxy"

n_subs_gal$panama_time_hms <- hms::hms(hours = n_subs_gal$panama_time)
n_subs_pres$group <- "Presidente"
n_subs_pres$panama_time_hms <- hms::hms(hours = n_subs_pres$panama_time)
n_subs_pres$sub5 <- NA

n_subs <- rbind(n_subs_gal, n_subs_pres)

#make plot for the number of subgroups for each hour of the day for both the groups

gg <- ggplot(n_subs, aes(panama_time_hms, n_groups, color = group))+
  geom_point(stat = "sum", position = position_dodge(width = 5000), alpha = 0.7) +
  # geom_violin(aes(color = group), trim = FALSE, position = position_dodge(0.5))+
  theme_classic()+ 
  scale_color_manual(values = c("Galaxy" = "darkolivegreen3", "Presidente" = "aquamarine3"))+
  labs(color = "Group") +
  #scale_x_binned(breaks = c(-Inf, 10, 14, Inf)) + #this only works if panama_time is not an hms object
  xlab("Time of day") +
  ylab("Number of subgroups") +
  NULL

gg

ggsave(filename = paste0(plot_dir, 'n_groups_time_of_day_dots.png'), plot = gg, width = 6, height = 4, dpi = 300)

