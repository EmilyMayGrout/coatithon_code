#This script is reading in the dataframes made in split_mechanics_exploration to plot the distance travelled by the larger and smaller subgroup for galaxy and Presidente combined for fissions and fusions seperately
library(ggplot2)
library(patchwork)

#load in the data


load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxy_fusion_grpsize.RData")
galaxy_fusion_df <- event_df
galaxy_fusion_df$Group <- "Galaxy"
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/presedente_fusion_grpsize.RData")
pres_fusion_df <- event_df
pres_fusion_df$Group <- "Presidente"
both_fusions <- rbind(galaxy_fusion_df, pres_fusion_df)



load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxy_fission_grpsize.RData")
galaxy_fission_df <- event_df
galaxy_fission_df$Group <- "Galaxy"
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/presedente_fission_grpsize.RData")
pres_fission_df <- event_df
pres_fission_df$Group <- "Presidente"
rm("event_df")
both_fissions <- rbind(galaxy_fission_df, pres_fission_df)


g1 <- ggplot(both_fusions, aes(x = Distance_Smaller_Group, 
                         y = Distance_Larger_Group, 
                         color = Group))+
  geom_point(alpha = 0.8, size = 2)+
  theme_classic()+
  labs(x = "Distance travelled by the smaller subgroup (m/s)", y = "Distance travelled by the larger subgroup (m/s)", title = "Fusions")+
  ylim(0, 80)+
  xlim(0, 80)+
  scale_color_manual(values=c("olivedrab2", "skyblue3"))+
  scale_alpha(guide = 'none')

g2 <- ggplot(both_fissions, aes(x = Distance_Smaller_Group, 
                          y = Distance_Larger_Group, 
                          color = Group))+
  geom_point(alpha = 0.8, size = 2, show.legend = F)+
  theme_classic()+
  labs(x = "Distance travelled by the smaller subgroup (m/s)", y = "Distance travelled by the larger subgroup (m/s)", title = "Fissions")+
  ylim(0, 80)+
  xlim(0, 80)+
  scale_color_manual(values=c( "olivedrab2", "skyblue3"))+
  scale_alpha(guide = 'none')

g2 + g1

ggsave(("C:/Users/egrout/Dropbox/coatithon/results/speeds_grpsize_fisfus.png"), width = 13, height = 5)

#how many events are smaller groups travelling more than bigger groups?
sum(both_fissions$Distance_Smaller_Group[both_fissions$Group == "Presidente"] > both_fissions$Distance_Larger_Group[both_fissions$Group == "Presidente"])
#37 events in Galaxy where smaller groups travel further
#11 for Presidente
sum(both_fissions$Distance_Smaller_Group[both_fissions$Group == "Presidente"] < both_fissions$Distance_Larger_Group[both_fissions$Group == "Presidente"])
#27 events in Galaxy where larger groups travel further
#25 for Presidente

sum(both_fusions$Distance_Smaller_Group[both_fusions$Group == "Galaxy"] > both_fusions$Distance_Larger_Group[both_fusions$Group == "Galaxy"])
#46 events in Galaxy where smaller groups travel to large group
#32 events in Presidente 
sum(both_fusions$Distance_Smaller_Group[both_fusions$Group == "Galaxy"] < both_fusions$Distance_Larger_Group[both_fusions$Group == "Galaxy"])
#13 in Galaxy where larger group travels further
#6 in Presidente



















