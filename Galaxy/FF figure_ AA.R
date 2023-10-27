# FF figure for Emily - this code was written by Alie Ashbury 
# 2023-10-24

library(tidyverse)
library(viridis)
library(plotly)

#read in the data
df <- read.csv("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/050122_gal_subgrouping.csv")
str(df)

#make the times into times, check it
df$ts <- as.POSIXct(df$ts)
str(df)
summary(df)

#prep the main subgrouping data for plotting
df_mod <- df %>%
  fill(sub.group, .direction = "downup") %>% #get rid of NAs
  mutate(new_sg = sub.group) %>% #make a new subgroup column in case we mess with it (unnecessary step)
  arrange(ts, coati_ids.name, sub.group) %>% #sort by time, coati id, and subgroup
  group_by(ts, sub.group) %>% #group by time and subgroup
  #make new columns...
  mutate(n_in_sg = n(), #number of coatis in each subgroup 
         seq_in_sg = seq_along(coati_ids.name), #number each coati in each subgroup
         sg_jitter = new_sg - (0.07/2*n_in_sg) + (seq_in_sg*0.07)) # this is a trixie but vital little pieces that
                                                                   # gives each coati in each subgroup a unique position,
                                                                   # the position is relative to the other coatis in the subroup,
                                                                   # and centered around the subgroup number


#plot it for the first time (then we will make some subgroup value adjustments)
p1 <- ggplot() +
  #plot a line for each individual, coloured by individual
  geom_line(data = df_mod,
               mapping = aes(x = ts,
                             y = sg_jitter,
                             group = as.factor(coati_ids.name),
                             color= as.factor(coati_ids.name)),
               alpha = 0.6,
               size = 0.7) +
  #plot a point for each individual, coloured by individual
  geom_point(data = df_mod,
             mapping = aes(x = ts,
                           y = sg_jitter,
                           color = as.factor(coati_ids.name)),
             size = 0.9,
             alpha = 0.9) +
  # use a color scale that is nicer than the default
  scale_color_manual(values = viridis(12)[1:11])+ 
  theme_classic() +
  # get rid of grid lines and extra lables etc that we don't need
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

p1 #look at it - notice where subgroup lines cross each other unnecessarily!
# if we add a subgroup position BELOW the main group, and pull some of the subgroup 2s down to there, it will be cleaner!
# use ggplotly to figure out exactly which points should be moved...
ggplotly(p1) #this allows you to hover over points with your cursor and see their aesthetic

# so now, we can use the above interactive plot to figure out which points need to be manually reset.
# we will move those points to a new subgroup that is below subgroup 1, ie subgroup 0 

df_mod2 <- df_mod #lets duplicated it before we mess with it, just in case we eff it up

#the following three little code chunks replace the df_mod2$new_sg for specific rows with 0
df_mod2$new_sg[which(df_mod2$coati_ids.name >=6 & 
                       df_mod2$coati_ids.name <= 9 &
                       df_mod2$ts == as.POSIXct("2022-01-05 12:20:00"))] <- 0

df_mod2$new_sg[which(df_mod2$coati_ids.name == 11 & 
                       (df_mod2$ts == as.POSIXct("2022-01-05 20:20:00") |
                          df_mod2$ts == as.POSIXct("2022-01-05 21:00:00") ))] <- 0

df_mod2$new_sg[which(df_mod2$sub.group != 1 & 
                       (df_mod2$ts >= as.POSIXct("2022-01-05 13:50:00") &
                          df_mod2$ts <= as.POSIXct("2022-01-05 15:40:00") ))] <- 0

#now we need to recalculate the special subgroup jitter values, since we've changed some of the subgroup values
# (same code as above basically)
df_mod2 <- df_mod2 %>%
  arrange(ts, coati_ids.name, sub.group) %>%
  group_by(ts, sub.group) %>%
  mutate(n_in_sg = n(),
         seq_in_sg = seq_along(coati_ids.name),
         sg_jitter = new_sg - (0.07/2*n_in_sg) + (seq_in_sg*0.07))

#let's make some little rectangles that can outline each subgroup
df_group_polys <- df_mod2 %>%
  group_by(ts, sub.group) %>% # group by timestamp and subgroup
  summarize(top_subgroup = max(sg_jitter) + 0.1, #take the highest jitter value and add a little buffer
            bottom_subgroup = min(sg_jitter) - 0.1, #take the lowest jitter value and add a little buffer
            left_subgroup = mean(ts) - 180, #take that timestamp and subtract a little buffer
            right_subgroup = mean(ts) + 180) # take that time stamp and add a little buffer
#so, above, if you want more white space between the subgroup outlines and the points, use bigger buffers

#let's make some little lines that will connect subgroup boxes at the same timestamp (otherwise it's hard to see how they line up)
df_ts_lines <- df_mod2 %>%
  group_by(ts) %>%
  summarize(top_subgroup= max(sg_jitter), #take the highest point at each timestamp
            bottom_subgroup = min(sg_jitter)) # take the lowest point at each timestamp

dev.off() #needed to reset plotting area after using ggplotly (above)

#ok now we plot it again
p2 <- ggplot() +
  #plot the little lines that link the subgroups
  geom_segment(data = df_ts_lines,
               mapping = aes(x = ts,
                             xend = ts,
                             y = bottom_subgroup,
                             yend = top_subgroup),
               color = "grey80",
               size = 0.3) +
  #plot little white rectangles over those lines (so that the lines don't connect points together across time stamps,
  # but instead only go from the edgest of the subgroup boxes)
  geom_rect(data = df_group_polys,
            mapping = aes(xmin = left_subgroup,
                          xmax = right_subgroup,
                          ymin = bottom_subgroup,
                          ymax = top_subgroup),
            fill = "white", #white fill to cover the lines plotted above
            color = "transparent") + #no outline
  #plot the lines for each coati
  geom_line(data = df_mod2, 
            mapping = aes(x = ts,
                          y = sg_jitter,
                          group = as.factor(coati_ids.name),
                          color = as.factor(coati_ids.name)),
            alpha = 0.3,
            size = 1) +
  #plot the points for each coati
  geom_point(data = df_mod2,
             mapping = aes(x = ts,
                           y = sg_jitter,
                           color = as.factor(coati_ids.name)),
             size = 0.6,
             alpha = 1) +
  #now plot the little subgroup rectangles again, but ON TOP of the lines & points, so that the outlines can
  # be clearly seen
  geom_rect(data = df_group_polys,
            mapping = aes(xmin = left_subgroup,
                          xmax = right_subgroup,
                          ymin = bottom_subgroup,
                          ymax = top_subgroup),
            fill = "transparent", #transparent fill so that they don't block otu the coati points and lines
            color = "grey40") + #dark grey outline
  #now choose some nice differentiated colors for each coati
  scale_color_manual(values = sample(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928'), #these are 11 distinct colors, they will be randomly assigned to the coatis
                                     11, replace = FALSE)) +
  theme_classic() +
  #get rid of all the extra stuff (note, this now is for a vertical plot. if you want a horizontal plot,
  # then replace the xs with ys and the ys with xs)
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
# +coord_flip() # flip it so that it's vertical

p2 #look at it

# save it
ggsave(filename = "FF line and dot subgroupings, vert.png",
       width = 2.5,
       height = 8,
       units = "in",
       dpi = 300,
       scale = 1.3)