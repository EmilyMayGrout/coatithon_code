#this script is for plotting the two groups that exhibit subgrouping behaviour the age/sex class of the individuals who are on their own

library(ggplot2)

#dataframes made in fission_fusion_galaxy_V1 and fission_fusion_presedente_V1
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/gal_alone_inds_level1.Rdata") 
load("C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/pres_alone_inds_level1.Rdata") 
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/trago_alone_inds_level1.Rdata") 

#add column for group name so can distinguish between the groups in the plot
coati_ids_alone_gal$group <- "Galaxy"
coati_ids_alone_pres$group <- "Presidente"
coati_ids_alone_trago$group <- "Trago"

#remove Tequila and Rum because they leave the group permenantly
coati_ids_alone_trago <- coati_ids_alone_trago[-c(4,6),]


gal_pres_trag_inds_alone <- rbind(coati_ids_alone_gal, coati_ids_alone_pres, coati_ids_alone_trago)

colors <- c("darkolivegreen3", "aquamarine3","hotpink3")

gg <- ggplot(aes(x = prop_alone, y = age_sex), data = gal_pres_trag_inds_alone) +
  xlab("Proportion of time alone (%)") +
  ylab("Age class") +
  geom_point(aes(color = factor(group)), shape = 16, size = 3, alpha = 0.7, position = position_jitter(width = 0.001, height = 0.25, seed = 312)) +
  facet_grid(vars(age_sex), scales = "free", space = "free") +
  scale_color_manual(values = colors) +
  theme_classic() +
  #guides(color = "none", shape = "none", alpha = "none") +  # Remove the legends
  theme(strip.text.y = element_blank())+# Remove the alpha legend 
  guides(color=guide_legend("Group:"))

gg

#for visualisation purposes, the data was jittered

ggsave(filename = 'C:/Users/egrout/Dropbox/coatithon/results/inds_alone_both_groups.png', plot = gg, width = 6, height = 6, dpi = 300)


