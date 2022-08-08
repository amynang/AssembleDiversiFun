library(tidyverse)
library(tidybayes)
library(see)
library(scico)
library(patchwork)
library(png)
library(magick)
library(cowplot)
set.seed(321)

properties = readRDS("properties_20220803_N_40_comp_40000cont.RData")

hh <- readPNG("C:/Users/aa21qeqa/Documents/AssembleDiversiFun/hh.png", 
               native = TRUE)
lh <- readPNG("C:/Users/aa21qeqa/Documents/AssembleDiversiFun/lh.png", 
               native = TRUE)
hi <- readPNG("C:/Users/aa21qeqa/Documents/AssembleDiversiFun/hi.png", 
               native = TRUE)
hl <- readPNG("C:/Users/aa21qeqa/Documents/AssembleDiversiFun/hl.png", 
               native = TRUE)
ll <- readPNG("C:/Users/aa21qeqa/Documents/AssembleDiversiFun/ll.png", 
               native = TRUE)

delta = properties %>% mutate_at(vars(plant.biomass), ~ifelse(.==0,NA,.)) %>% 
  group_by(ID) %>%
  summarise(plant.rich = first(plant.rich),
            deltaAnimals = plant.biomass[scenario == 'low for animals - high inter for plants'] - plant.biomass[scenario == 'high for animals - high inter for plants'],
            deltaPlants1 = plant.biomass[scenario == 'high for animals - high intra for plants'] - plant.biomass[scenario == 'high for animals - high inter for plants'],
            deltaPlants2 = plant.biomass[scenario == 'high for animals - low for plants'] - plant.biomass[scenario == 'high for animals - high inter for plants'],
            deltaBoth = plant.biomass[scenario == 'low for animals - low for plants'] - plant.biomass[scenario == 'high for animals - high inter for plants'])

a = ggplot(delta, aes(log2(plant.rich), deltaAnimals)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  #theme_modern()+
  theme_modern() #+
  #inset_element(p = img,
  #              left = 0.5,
  #              bottom = 0.5,
  #              right = 0.9,
  #              top = 0.9)

 # a = ggdraw() +
 #  draw_plot(a) +
 #  draw_image(img, x = .35, y = -.2, scale = .2)
  


p1 = ggplot(delta, aes(log2(plant.rich), deltaPlants1)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  #theme_modern()+
  theme_modern()

p2 = ggplot(delta, aes(log2(plant.rich), deltaPlants2)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  #theme_modern()+
  theme_modern()

b = ggplot(delta, aes(log2(plant.rich), deltaBoth)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .05) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  theme_modern()


baseline =  ggplot(properties[which(properties$scenario == "high for animals - high inter for plants"
                                    & properties$plant.end!=0),], 
                   aes(log2(plant.rich), plant.biomass)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .1) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  theme_modern()

#(baseline + (a+p1)/(p2+b)) + plot_annotation(tag_levels = 'a')


# left, baseline
cow_a = plot_grid(baseline,
                  ncol = 1,
                  labels = "a",
                  label_size = 15)
# right, arrange the differences to baseline
cow_bcde = plot_grid(a, p1, p2, b,
                     ncol = 2,
                     labels = c("b","c","d","e"),
                     label_size = 15)
# put the two together
cow_abcde = plot_grid(cow_a, cow_bcde,
                      ncol = 2)
# add the inset scenario pictures
cow_final <- ggdraw() +
  draw_plot(cow_abcde) +
  draw_image(image = hh, 
             x = -.1, 
             y =  -.28, 
             scale = .2) +
  draw_image(image = lh, 
             x = .18, 
             y =  .18, 
             scale = 0.15) +
  draw_image(image = hl, 
             x = .18, 
             y =  -.32, 
             scale = .15) +
  draw_image(image = hi, 
             x = .453, 
             y =  .18, 
             scale = .15) +
  draw_image(image = ll, 
             x = .43, 
             y =  -.32, 
             scale = .15)
# save the plot
ggsave("plantbiomass3.png",
       scale = 2,
       width = 173,
       height = 173/2,
       units = "mm",
       dpi = 600)
########################## Primary productivity ################################

delta = properties %>% mutate_at(vars(primary.productivity), ~ifelse(.==0,NA,.)) %>% 
  group_by(ID) %>%
  summarise(plant.rich = first(plant.rich),
            deltaAnimals = primary.productivity[scenario == 'low for animals - high inter for plants'] - primary.productivity[scenario == 'high for animals - high inter for plants'],
            deltaPlants1 = primary.productivity[scenario == 'high for animals - high intra for plants'] - primary.productivity[scenario == 'high for animals - high inter for plants'],
            deltaPlants2 = primary.productivity[scenario == 'high for animals - low for plants'] - primary.productivity[scenario == 'high for animals - high inter for plants'],
            deltaBoth = primary.productivity[scenario == 'low for animals - low for plants'] - primary.productivity[scenario == 'high for animals - high inter for plants'],
            deltaBothVSPlants2 = primary.productivity[scenario == 'low for animals - low for plants'] - primary.productivity[scenario == 'high for animals - low for plants'])


a = ggplot(delta, aes(log2(plant.rich), deltaAnimals)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-productivity"),
       x = "") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c('2', '4', '8', '16')) +
  theme_modern()


p1 = ggplot(delta, aes(log2(plant.rich), deltaPlants1)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-productivity"),
       x = "") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c('2', '4', '8', '16')) +
  theme_modern()

p2 = ggplot(delta, aes(log2(plant.rich), deltaPlants2)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-productivity"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()

b = ggplot(delta, aes(log2(plant.rich), deltaBoth)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-productivity"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()


baseline =  ggplot(properties[which(properties$scenario == "high for animals - high inter for plants"
                                    & properties$plant.end!=0),], 
                   aes(log2(plant.rich), primary.productivity)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = "Primary productivity",
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()

#(baseline + (a+p1)/(p2+b))  + plot_annotation(tag_levels = 'a')

# left, baseline
cow_a = plot_grid(baseline,
                  ncol = 1,
                  labels = "a",
                  label_size = 15)
# right, arrange the differences to baseline
cow_bcde = plot_grid(a, p1, p2, b,
                     ncol = 2,
                     labels = c("b","c","d","e"),
                     label_size = 15)
# put the two together
cow_abcde = plot_grid(cow_a, cow_bcde,
                      ncol = 2)
# add the inset scenario pictures
cow_final <- ggdraw() +
  draw_plot(cow_abcde) +
  draw_image(image = hh, 
             x = -.1, 
             y =  -.28, 
             scale = .2) +
  draw_image(image = lh, 
             x = .18, 
             y =  .18, 
             scale = 0.15) +
  draw_image(image = hl, 
             x = .18, 
             y =  -.32, 
             scale = .15) +
  draw_image(image = hi, 
             x = .453, 
             y =  .18, 
             scale = .15) +
  draw_image(image = ll, 
             x = .43, 
             y =  -.32, 
             scale = .15)
# save the plot
ggsave("primaryprod2.png",
       scale = 2,
       width = 180,
       height = 180/2,
       units = "mm",
       dpi = 600)

ggplot(delta, aes(log2(plant.rich), deltaBothVSPlants2)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-productivity"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c('2', '4', '8', '16')) +
  theme_modern()

ggsave("animalsVSjustplants.png",
       scale = 2,
       width = 180/2,
       height = 180/2,
       units = "mm",
       dpi = 600)

############################## Herbivory pressure ##############################

delta = properties %>% mutate_at(vars(herbivory.pressure), ~ifelse(.==0,NA,.)) %>% 
  group_by(ID) %>%
  summarise(plant.rich = first(plant.rich),
            deltaAnimals = herbivory.pressure[scenario == 'low for animals - high inter for plants'] - herbivory.pressure[scenario == 'high for animals - high inter for plants'],
            deltaPlants1 = herbivory.pressure[scenario == 'high for animals - high intra for plants'] - herbivory.pressure[scenario == 'high for animals - high inter for plants'],
            deltaPlants2 = herbivory.pressure[scenario == 'high for animals - low for plants'] - herbivory.pressure[scenario == 'high for animals - high inter for plants'],
            deltaBoth = herbivory.pressure[scenario == 'low for animals - low for plants'] - herbivory.pressure[scenario == 'high for animals - high inter for plants'])

a = ggplot(delta, aes(log2(plant.rich), deltaAnimals)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-pressure"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()


p1 = ggplot(delta, aes(log2(plant.rich), deltaPlants1)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-pressure"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()

p2 = ggplot(delta, aes(log2(plant.rich), deltaPlants2)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-pressure"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()

b = ggplot(delta, aes(log2(plant.rich), deltaBoth)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-pressure"),
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()


baseline =  ggplot(properties[which(properties$scenario == "high for animals - high inter for plants"
                                    & properties$plant.end!=0),], 
                   aes(log2(plant.rich), herbivory.pressure)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = "Herbivory pressure on plants",
       x = "Initial plant richness") +
    scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                       labels = c('2', '4', '8', '16')) +   
  theme_modern()

#(baseline + (a+p1)/(p2+b)) + plot_annotation(tag_levels = 'a')

# left, baseline
cow_a = plot_grid(baseline,
                  ncol = 1,
                  labels = "a",
                  label_size = 15)
# right, arrange the differences to baseline
cow_bcde = plot_grid(a, p1, p2, b,
                     ncol = 2,
                     labels = c("b","c","d","e"),
                     label_size = 15)
# put the two together
cow_abcde = plot_grid(cow_a, cow_bcde,
                      ncol = 2)
# add the inset scenario pictures
cow_final <- ggdraw() +
  draw_plot(cow_abcde) +
  draw_image(image = hh, 
             x = -.1, 
             y =  .38, 
             scale = .2) +
  draw_image(image = lh, 
             x = .18, 
             y =  .18, 
             scale = 0.15) +
  draw_image(image = hl, 
             x = .18, 
             y =  -.32, 
             scale = .15) +
  draw_image(image = hi, 
             x = .453, 
             y =  .18, 
             scale = .15) +
  draw_image(image = ll, 
             x = .43, 
             y =  -.32, 
             scale = .15)
# save the plot
ggsave("herbivory.png",
       scale = 2,
       width = 180,
       height = 180/2,
       units = "mm",
       dpi = 600)

####################### Control of herbivores by predators ######################

delta = properties %>% #mutate_at(vars(herbivory.control), ~ifelse(.==0,NA,.)) %>% 
  filter(herbivory.control < quantile(herbivory.control, probs = .99, na.rm = T)) %>% 
  group_by(ID) %>%
  summarise(plant.rich = first(plant.rich),
            deltaAnimals = herbivory.control[scenario == 'low for animals - high inter for plants'] - herbivory.control[scenario == 'high for animals - high inter for plants'],
            deltaPlants1 = herbivory.control[scenario == 'high for animals - high intra for plants'] - herbivory.control[scenario == 'high for animals - high inter for plants'],
            deltaPlants2 = herbivory.control[scenario == 'high for animals - low for plants'] - herbivory.control[scenario == 'high for animals - high inter for plants'],
            deltaBoth = herbivory.control[scenario == 'low for animals - low for plants'] - herbivory.control[scenario == 'high for animals - high inter for plants'],
            deltaBothVSPlants2 = primary.productivity[scenario == 'low for animals - low for plants'] - primary.productivity[scenario == 'high for animals - low for plants'])

a = ggplot(delta, aes(log2(plant.rich), deltaAnimals)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                     labels = c('2', '4', '8', '16')) + 
  theme_modern()


p1 = ggplot(delta, aes(log2(plant.rich), deltaPlants1)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                     labels = c('2', '4', '8', '16')) + 
  theme_modern()

p2 = ggplot(delta, aes(log2(plant.rich), deltaPlants2)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                     labels = c('2', '4', '8', '16')) + 
  theme_modern()

b = ggplot(delta, aes(log2(plant.rich), deltaBoth)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                     labels = c('2', '4', '8', '16')) + 
  theme_modern()


baseline =  ggplot(properties[which(properties$scenario == "high for animals - high inter for plants"
                                    & properties$plant.end!=0 &
                                      properties$herbivory.control < quantile(properties$herbivory.control, probs = .99, na.rm = T)),], 
                   aes(log2(plant.rich), herbivory.control)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = "Control of herivores by predators",
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),                      
                     labels = c('2', '4', '8', '16')) + 
  theme_modern()

#(baseline + (a+p1)/(p2+b)) + plot_annotation(tag_levels = 'a')

# left, baseline
cow_a = plot_grid(baseline,
                  ncol = 1,
                  labels = "a",
                  label_size = 15)
# right, arrange the differences to baseline
cow_bcde = plot_grid(a, p1, p2, b,
                     ncol = 2,
                     labels = c("b","c","d","e"),
                     label_size = 15)
# put the two together
cow_abcde = plot_grid(cow_a, cow_bcde,
                      ncol = 2)
# add the inset scenario pictures
cow_final <- ggdraw() +
  draw_plot(cow_abcde) +
  draw_image(image = hh, 
             x = -.1, 
             y =  .38, 
             scale = .2) +
  draw_image(image = lh, 
             x = .18, 
             y =  .18, 
             scale = 0.15) +
  draw_image(image = hl, 
             x = .18, 
             y =  -.32, 
             scale = .15) +
  draw_image(image = hi, 
             x = .453, 
             y =  .18, 
             scale = .15) +
  draw_image(image = ll, 
             x = .43, 
             y =  -.32, 
             scale = .15)
# save the plot
ggsave("control.png",
       scale = 2,
       width = 180,
       height = 180/2,
       units = "mm",
       dpi = 600)

ggplot(delta, aes(log2(plant.rich), deltaBothVSPlants2)) +
  geom_point(position = position_jitter(width = .3),
             color = "grey",
             alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(1,3),
                     point_size = 3
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(.5,2),
                     point_size = 2
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  scale_x_continuous(breaks = c(1, 2, 3, 4),
                     labels = c('2', '4', '8', '16')) +
  theme_modern()

ggsave("animalsVSjustplants-control.png",
       scale = 2,
       width = 180/2,
       height = 180/2,
       units = "mm",
       dpi = 600)

ggplot(properties[which(properties$scenario == "high for animals - high inter for plants"
                        & properties$plant.end!=0 &
                          properties$herbivory.control < quantile(properties$herbivory.control, probs = .99, na.rm = T)),], 
       aes(log2(plant.rich), herbivory.control)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .1) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  theme_modern()

ggplot(properties[which(properties$scenario == "low for animals - low for plants"
                        & properties$plant.end!=0 &
                          properties$herbivory.control < quantile(properties$herbivory.control, probs = .99, na.rm = T)),], 
       aes(log2(plant.rich), herbivory.control)) +
  geom_point(position = position_jitter(width = .3),
             alpha = .1) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "white",
                     interval_size_range = c(2,4),
                     point_size = 6
  ) +
  stat_pointinterval(point_interval = "mean_qi",
                     .width = c(0.75, .95), color = "black",
                     interval_size_range = c(1,3),
                     point_size = 5
  ) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#b6543e") +
  labs(y = "Herbivory pressure on plants",
       x = "Initial plant richness") +
  theme_classic()

############################## Species persistence #############################
ggplot(properties, 
       aes(x = log2(plant.rich), 
           y=plant.end/plant.rich, 
           color = scenario)) +
  geom_point(alpha = .02,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting plant species",
       x = "Initial plant richness")

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=animals.end/40, 
           color = scenario)) +
  geom_point(alpha = .02,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting animal species",
       x = "Initial plant richness")



ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=animal.biomass, 
           color = scenario)) +
  geom_point(alpha = .5,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting animal species",
       x = "Initial plant richness")



ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=animal.biomass, 
           color = scenario)) +
  geom_point(alpha = .5,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting animal species",
       x = "Initial plant richness")

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=primary.productivity, 
           color = scenario)) +
  geom_point(alpha = .5,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting animal species",
       x = "Initial plant richness")

ggplot(properties[which(properties$plant.end!=0 &
                          properties$herbivory.control < quantile(properties$herbivory.control, probs = .99, na.rm = T)),], 
       aes(x = log2(plant.rich), 
           y=(herbivory.control), 
           color = scenario)) +
  geom_point(alpha = .5,
             position = position_jitterdodge(dodge.width = .9)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), 
             linetype = "dashed",
             color = "grey") +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               aes(fill = scenario, 
                   shape = scenario),
               color = "black",
               position = position_dodge(width = .9)) +
  scale_shape_manual(values=c(21, 21, 21,
                              22, 22, 22)) +
  scale_fill_manual(values = c("#b6543e","#dda64a","#618d7a",
                               "#b6543e","#dda64a","#618d7a")) +
  scale_color_manual(values = c("#b6543e","#dda64a","#618d7a",
                                "#b6543e","#dda64a","#618d7a")) +
  theme_modern() + 
  labs(y = "Proportion of persisting animal species",
       x = "Initial plant richness")

tryme = properties %>% group_by(plant.rich, scenario) %>% 
  summarize(lost = (sum(plant.end == 0)/1000)*100)



props = properties %>% mutate(prop.plant = plant.biomass/(animal.biomass+plant.biomass)
                              ,
                              prop.herb  = herb.biomass /(animal.biomass+plant.biomass)
                              ,
                              prop.promn = (pred.biomass+omn.biomass)/(animal.biomass+plant.biomass)
                              ) %>% 
  select(c(1:2,19:21)) %>% pivot_longer(3:5, 
                                        names_to = "level",
                                        values_to = "prop.mass") %>% 
  ggplot(aes(log2(plant.rich), prop.mass, color = level)) +
  #geom_point(position = position_jitter(.3),
  #           alpha = .1) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_color_manual(values = c("#b6543e","#618d7a","#dda64a")) +
  facet_wrap(~scenario) +
  theme_classic()

properties %>% mutate(prop.omn = omn.biomass/(animal.biomass)
                      ,
                      prop.herb  = herb.biomass /(animal.biomass)
                      ,
                      prop.pred = pred.biomass/(animal.biomass)
) %>% 
  select(c(1:2,19:21)) %>% pivot_longer(3:5, 
                                        names_to = "level",
                                        values_to = "prop.mass") %>% 
  ggplot(aes(log2(plant.rich), prop.mass, color = level)) +
  geom_point(position = position_jitter(.3),
             alpha = .05) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_color_manual(values = c("#b6543e","#618d7a","#dda64a")) +
  facet_wrap(~scenario) +
  theme_classic()




