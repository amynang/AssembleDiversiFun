library(tidyverse)
library(tidybayes)
library(see)
library(scico)
#library(patchwork)
library(png)
library(magick)
library(cowplot)
library(grid)
library(ggthemes)
set.seed(321)

#properties = readRDS("properties_20220810_N_40_comp_40000qfixed.RData")
properties = readRDS("properties_2-16.RData")
properties$animal.competition = rep(rep(c("high","low"), each = 3), 4500)
properties$plant.competition = rep(rep(c("high inter",
                                         "high intra",
                                         "low"), 2), 4500)

s1 = "high for animals - high inter for plants"
s2 = "high for animals - high intra for plants"
s3 = "high for animals - low for plants"
s4 = "low for animals - high inter for plants"
s5 = "low for animals - high intra for plants"
s6 = "low for animals - low for plants"

leg = readPNG("./legend.png", 
              native = TRUE)
i2_1 = readPNG("./2_1.png", 
               native = TRUE)
i3_1 = readPNG("./3_1.png", 
               native = TRUE)
i6_3 = readPNG("./6_3.png", 
               native = TRUE)
i4_1 = readPNG("./4_1.png", 
               native = TRUE)
i5_2 = readPNG("./5_2.png", 
               native = TRUE)


delta = properties %>% 
  group_by(ID) %>%
  summarise(plant.rich = first(plant.rich),
            d2_1 = herbivory.control[scenario == s2] - herbivory.control[scenario == s1],
            d3_1 = herbivory.control[scenario == s3] - herbivory.control[scenario == s1],
            d4_1 = herbivory.control[scenario == s4] - herbivory.control[scenario == s1],
            d5_1 = herbivory.control[scenario == s5] - herbivory.control[scenario == s1],
            d6_1 = herbivory.control[scenario == s6] - herbivory.control[scenario == s1],
            d3_2 = herbivory.control[scenario == s3] - herbivory.control[scenario == s2],
            d4_2 = herbivory.control[scenario == s4] - herbivory.control[scenario == s2],
            d5_2 = herbivory.control[scenario == s5] - herbivory.control[scenario == s2],
            d6_2 = herbivory.control[scenario == s6] - herbivory.control[scenario == s2],
            d4_3 = herbivory.control[scenario == s4] - herbivory.control[scenario == s3],
            d5_3 = herbivory.control[scenario == s5] - herbivory.control[scenario == s3],
            d6_3 = herbivory.control[scenario == s6] - herbivory.control[scenario == s3],
            d5_4 = herbivory.control[scenario == s5] - herbivory.control[scenario == s4],
            d6_4 = herbivory.control[scenario == s6] - herbivory.control[scenario == s4],
            d6_5 = herbivory.control[scenario == s6] - herbivory.control[scenario == s5])

prop = delta %>% group_by(plant.rich) %>% 
  summarise(pp2_1 = sum(d2_1>0, na.rm = T)/300,
            pp3_1 = sum(d3_1>0, na.rm = T)/300,
            pp4_1 = sum(d4_1>0, na.rm = T)/300,
            pp5_1 = sum(d5_1>0, na.rm = T)/300,
            pp6_1 = sum(d6_1>0, na.rm = T)/300,
            pp3_2 = sum(d3_2>0, na.rm = T)/300,
            pp4_2 = sum(d4_2>0, na.rm = T)/300,
            pp5_2 = sum(d5_2>0, na.rm = T)/300,
            pp6_2 = sum(d6_2>0, na.rm = T)/300,
            pp4_3 = sum(d4_3>0, na.rm = T)/300,
            pp5_3 = sum(d5_3>0, na.rm = T)/300,
            pp6_3 = sum(d6_3>0, na.rm = T)/300,
            pp5_4 = sum(d5_4>0, na.rm = T)/300,
            pp6_4 = sum(d6_4>0, na.rm = T)/300,
            pp6_5 = sum(d6_5>0, na.rm = T)/300) %>% 
  pivot_longer(cols = 2:16,
               names_to = "difference",
               values_to = "p")







a = ggplot(properties, 
           aes((plant.rich), herbivory.control,
               color = plant.competition)) +
  geom_smooth(aes(linetype = animal.competition), se = T,
              method = "gam", formula = y ~ s(x, bs = "cs")
  ) +
  labs(y = "Herbivore control",
       x = "Initial plant richness") +
  guides(linetype = guide_legend(override.aes= list(color = "black"))) +
  #coord_equal(ratio = 16/3.5) +
  scale_color_manual(values = c("#dda64a","#b6543e","#618d7a")) +
  scale_linetype_manual(values = c("solid","twodash")) +
  theme_modern() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2)),
        plot.margin = unit(c(5.5, 5.5, 5.5, 15.5), "pt"),
        legend.position = "none")


p = ggplot(delta[,1:4] %>% pivot_longer(3:4,
                                        names_to = "which",
                                        values_to = "value"), 
           aes((plant.rich), value)) +
  geom_point(aes(color = which),
             size = 1,
             position = position_jitter(width = .3),
             #color = "grey",
             alpha = .25,
             show.legend = F) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "white",
  #                   interval_size_range = c(1,3),
  #                   point_size = 3
  #) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "black",
  #                   interval_size_range = c(.5,2),
  #                   point_size = 2
  #) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#dda64a", size = 1.25) +
  geom_smooth(#color = "#618d7a"
    aes(color = which), size = 1.25,
    show.legend = F) +
  scale_color_manual(values = c("#b6543e","#618d7a")) +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  coord_equal(ratio = 16/.8) +
  theme_modern() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2), margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")


a1 = ggplot(delta[,c(1:2,5)], 
            aes((plant.rich), d4_1)) +
  geom_point(#aes(color = which),
    size = 1,
    position = position_jitter(width = .3),
    color = "grey",
    alpha = .25) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "white",
  #                   interval_size_range = c(1,3),
  #                   point_size = 3
  #) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "black",
  #                   interval_size_range = c(.5,2),
  #                   point_size = 2
  #) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#dda64a", linewidth = 1) +
  geom_smooth(color = "#dda64a", linetype = "twodash", linewidth = 1) +
  #facet_wrap(~which, ncol = 1) +
  #scale_color_manual(values = c("#b6543e","#618d7a")) +
  labs(y = expression(delta~"-control"),
       x = " ") +
  ylim(-1,1) +
  coord_equal(ratio = 16/2) +
  theme_modern() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2), margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")

a2 = ggplot(delta[,c(1:2,10)], 
            aes((plant.rich), d5_2)) +
  geom_point(#aes(color = which),
    size = 1,
    position = position_jitter(width = .3),
    color = "grey",
    alpha = .25) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "white",
  #                   interval_size_range = c(1,3),
  #                   point_size = 3
  #) +
  #stat_pointinterval(point_interval = "mean_qi",
  #                   .width = c(0.75, .95), color = "black",
  #                   interval_size_range = c(.5,2),
  #                   point_size = 2
  #) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#b6543e", linewidth = 1) +
  geom_smooth(color = "#b6543e", linetype = "twodash", linewidth = 1) +
  #facet_wrap(~which, ncol = 1) +
  #scale_color_manual(values = c("#b6543e","#618d7a")) +
  labs(y = expression(delta~"-control"),
       x = "Initial plant richness") +
  ylim(-1,1) +
  coord_equal(ratio = 16/2) +
  theme_modern() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2), margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")

a3 = ggplot(delta[,c(1:2,14)], 
            aes((plant.rich), d6_3)) +
  geom_point(#aes(color = which),
    size = 1,
    position = position_jitter(width = .3),
    color = "grey",
    alpha = .25) +
  geom_hline(yintercept = 0, linetype = "solid", color = "#618d7a", linewidth = 1) +
  geom_smooth(color = "#618d7a", linetype = "twodash", linewidth = 1) +
  labs(y = expression(delta~"-control"),
       x = " ") +
  ylim(-1,1) +
  coord_equal(ratio = 16/2) +
  theme_modern() +
  theme(axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2), margin = margin(t = 0, r = 0, b = 0, l = 0)),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")



cow_p = plot_grid(NULL, p, NULL,
                  ncol = 1,
                  labels = c("","b",""),
                  label_size = 15,
                  rel_heights = c(1/4,2/4,1/4))
cow_a = plot_grid(a3, a1, a2,
                  ncol = 1,
                  labels = c("c","",""),
                  label_size = 15)

cow_prod = plot_grid(cow_p, cow_a,
                     ncol = 2,
                     labels = c("",""),
                     label_size = 15,
                     rel_widths = c(1,.75))

cow_c = plot_grid(NULL, a, NULL,
                  ncol = 1,
                  labels = c("","a",""),
                  label_size = 15,
                  rel_heights = c(1/4,2/4,1/4))
cow_pres = plot_grid(cow_c, cow_p, NULL, cow_a, NULL,
                     ncol = 5,
                     labels = c("","","",""),
                     label_size = 15,
                     rel_widths = c(1, 1, .1, .75, .1))


################################ Inset plots ###################################

ggplot(prop,aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey") + 
  facet_wrap(~difference, scales='free') +
  theme_tufte() +
  theme(axis.line=element_line()) + 
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

p2_1 = ggplot(prop[prop$difference == "pp2_1",],aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey", se = F) + 
  theme_tufte() +
  theme(axis.line=element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

p3_1 = ggplot(prop[prop$difference == "pp3_1",],aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey", se = F) + 
  #facet_wrap(~difference, scales='free') +
  theme_tufte() +
  theme(axis.line=element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

p4_1 = ggplot(prop[prop$difference == "pp4_1",],aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey", se = F) + 
  #facet_wrap(~difference, scales='free') +
  theme_tufte() +
  theme(axis.line=element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

p5_2 = ggplot(prop[prop$difference == "pp5_2",],aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey", se = F) + 
  #facet_wrap(~difference, scales='free') +
  theme_tufte() +
  theme(axis.line=element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

p6_3 = ggplot(prop[prop$difference == "pp6_3",],aes(plant.rich, p)) +
  geom_point() +
  geom_hline(yintercept = c(.25,.5,.75), color = "grey", linetype="dashed") +
  geom_smooth(color = "darkgrey", se = F) + 
  #facet_wrap(~difference, scales='free') +
  theme_tufte() +
  theme(axis.line=element_line(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_x_continuous(limits=c(1,16)) + 
  scale_y_continuous(limits=c(0,1))

ggdraw(cow_pres + 
         draw_plot(p2_1, 
                   x = .075, y = -.255,
                   width = .75, height = .75, 
                   scale = .23) +
         draw_plot(p3_1, 
                   x = .245, y = -.255,
                   width = .75, height = .75, 
                   scale = .23)+
         draw_plot(p6_3, 
                   x = .52, y = .428,
                   width = .75, height = .75, 
                   scale = .155)+
         draw_plot(p4_1, 
                   x = .52, y = .0935,
                   width = .75, height = .75, 
                   scale = .155)+
         draw_plot(p5_2, 
                   x = .52, y = -.24,
                   width = .75, height = .75, 
                   scale = .155) +
         draw_image(leg,
                    x = -.4, y = -.4,
                    scale = .2) +
         draw_image(i3_1,
                    x = .2, y = .15,
                    scale = .08) +
         draw_image(i2_1,
                    x = .196, y = .02,
                    scale = .075) +
         draw_image(i6_3,
                    x = .467, y = .4,
                    scale = .055) +
         draw_image(i4_1,
                    x = .467, y = .065,
                    scale = .055) +
         draw_image(i5_2,
                    x = .467, y = -.265,
                    scale = .055) +
         draw_grob(rectGrob(x = .45, y = .122,
                            width = .17,
                            height = .17,
                            gp = gpar(col = "#b6543e",
                                      fill = alpha("white", 0),
                                      alpha = 1,
                                      lwd = 3))) +
         draw_grob(rectGrob(x = .622, y = .122,
                            width = .17,
                            height = .17,
                            gp = gpar(col = "#618d7a",
                                      fill = alpha("white", 0),
                                      alpha = 1,
                                      lwd = 3))) +
         annotate(geom = "curve", # to red
                  size = .5,
                  x = .495, xend = .445,
                  y = .495, yend = .208,
                  curvature = 0) +
         annotate(geom = "curve", # to green
                  size = .5,
                  x = .595, xend = .625,
                  y = .575, yend = .208,
                  curvature = 0)
) 


ggsave("Figure3.png",
       bg = "white",
       scale = 1.9,
       width = 180,
       height = 180*(2/3),
       units = "mm",
       dpi = 300)
