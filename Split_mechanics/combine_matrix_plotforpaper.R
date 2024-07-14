#combining the galaxy and presidente matrices for plot in paper

load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxy_matrix.RData") #saved in split_mechanics_exploration
contingency_df_gal <- contingency_df
contingency_df_gal$group <- "Galaxy"
contingency_df_gal <- contingency_df_gal[contingency_df_gal$Var2 == "many/many",]
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/presedente_matrix.RData")
contingency_df_pres <- contingency_df
contingency_df_pres$group <- "Presidente"
contingency_df_pres <- contingency_df_pres[contingency_df_pres$Var2 == "many/many",]

contingency_df <- rbind(contingency_df_gal, contingency_df_pres)

all <- ggplot(contingency_df, aes(y = Var1, x = as.factor(group), fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "royalblue3") +
  labs(x = "Subgroup",
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

ggsave("C:/Users/egrout/Dropbox/coatithon/results/fission_matrix_bothgrps.png", all, width = 8, height = 8)




#make same matrix for fusions:


load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxyfusion_matrix.RData") #saved in split_mechanics_exploration
contingency_df_gal <- contingency_df
contingency_df_gal$group <- "Galaxy"
contingency_df_gal <- contingency_df_gal[contingency_df_gal$Var2 == "many/many",]
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/presedentefusion_matrix.RData")
contingency_df_pres <- contingency_df
contingency_df_pres$group <- "Presidente"
contingency_df_pres <- contingency_df_pres[contingency_df_pres$Var2 == "many/many",]

contingency_df <- rbind(contingency_df_gal, contingency_df_pres)


fus <- ggplot(contingency_df, aes(y = Var1, x = as.factor(group), fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "darkorchid3") +
  labs(x = "Subgroup",
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

fus

ggsave("C:/Users/egrout/Dropbox/coatithon/results/fusion_matrix_bothgrps.png", fus, width = 8, height = 8)














