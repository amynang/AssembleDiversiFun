library(tidyverse)


results = readRDS("results_20220721_N_40_comp_40000.RData")


# soll = as.data.frame(results[[36]][[1]])
# # colnames(soll) = c("time",
# #                    # paste0("nut_",1:2),
# #                    colnames(reg.loc[[n]][[s]]))
# 
# solll = soll %>% pivot_longer(!time, names_to = "species", values_to = "biomass")
# solll$taxon = substr(solll$species, 6, 8)
# 
# # plot results
# ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
#   geom_point() 




# w = 1
# outflux = results[[w]][[2]]$F#[i,j]
# outflux[,] = NA
# for (i in 1:nrow(results[[w]][[2]]$F)) {
#   for (j in 1:ncol(results[[w]][[2]]$F)) {
#     
#     outflux[i,j] = (results[[w]][[2]]$X[2+j]* #mass specific metabolic rate
#                     results[[w]][[2]]$max_feed[j]* # maximum feeding rate rel. to met. rate
#                     results[[w]][[1]][201,3+j]* # biomass at last time-step
#                     results[[w]][[2]]$F[i,j])/ # functional responses at last timestep
#                     results[[w]][[2]]$e[i] # assimilation efficiency
#   }
# }
# colnames(outflux) = colnames(results[[w]][[1]])[-(1:3)]
# rownames(outflux) = colnames(results[[w]][[1]])[-1]
# which(colSums(outflux)==0)
# sum(outflux[1:2,],na.rm = T)
# 
# 
# plant.growth = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
# for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
#   for (j in 1:dim(results[[w]][[2]]$alpha)[1]) {
#     plant.growth[i] = 1 - (results[[w]][[2]]$alpha[i,j]*
#                            results[[w]][[1]][201,1+j]/
#                            results[[w]][[2]]$K)
#   }
# }
# plant.prod = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
# for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
#   plant.prod[i] = results[[w]][[2]]$r[i]*
#                   plant.growth[i]*
#                   results[[w]][[1]][201,1+i]
#   
# }
# 
# 
# 
# w=1
properties = vector(mode = "list", length = length(results))

t.end = 201

for (w in 21619:length(results)) {
  
  # turn biomass of extinct species to 0
  results[[w]][[1]][t.end,] = ifelse(results[[w]][[1]][t.end,]<=1e-6, 
                                     0, 
                                     results[[w]][[1]][t.end,])
  properties[[w]]$scenario = NA
  
  properties[[w]]$plant.rich = length(grep("plant", colnames(results[[w]][[1]])))
  
  properties[[w]]$plant.end = length(which(results[[w]][[1]][t.end, grep("plant", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$animals.end = length(which(results[[w]][[1]][t.end, grep("herb|omn|pred", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$herb.end = length(which(results[[w]][[1]][t.end, grep("herb", colnames(results[[w]][[1]]))]>0))
  properties[[w]]$omn.end  = length(which(results[[w]][[1]][t.end, grep("omn",  colnames(results[[w]][[1]]))]>0))
  properties[[w]]$pred.end = length(which(results[[w]][[1]][t.end, grep("pred", colnames(results[[w]][[1]]))]>0))
  
  
  # mean number of predatory links of herbivores (incl. to omnivores)
  properties[[w]]$mean.predofherb = mean(rowSums(vegan::decostand(results[[w]][[2]]$F[which(results[[w]][[1]][t.end, 
                                                                                                              grep("herb", 
                                                                                                                   colnames(results[[w]][[1]]))]>0),], 
                                                                  "pa")))
  
  properties[[w]]$plant.biomass  = sum(results[[w]][[1]][t.end, grep("plant",         colnames(results[[w]][[1]]))])
  properties[[w]]$animal.biomass = sum(results[[w]][[1]][t.end, grep("herb|omn|pred", colnames(results[[w]][[1]]))])
  properties[[w]]$herb.biomass   = sum(results[[w]][[1]][t.end, grep("herb",          colnames(results[[w]][[1]]))])
  properties[[w]]$omn.biomass    = sum(results[[w]][[1]][t.end, grep("omn",           colnames(results[[w]][[1]]))])
  properties[[w]]$pred.biomass   = sum(results[[w]][[1]][t.end, grep("pred",          colnames(results[[w]][[1]]))])
  
  # calculate growth rate for each plant species
  plant.growth = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
  for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
    for (j in 1:dim(results[[w]][[2]]$alpha)[1]) {
      plant.growth[i] = 1 - (results[[w]][[2]]$alpha[i,j]* # competition from species j
                             results[[w]][[1]][t.end,1+j]/   # biomass of species j
                             results[[w]][[2]]$K)          # carrying capacity of species i
    }
  }
  #calculate plant productivity for each plant species
  plant.prod = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
  for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
    plant.prod[i] = results[[w]][[2]]$r[i]*
                    plant.growth[i]*
                    results[[w]][[1]][t.end,1+i]
  }
  
  properties[[w]]$primary.productivity = sum(plant.prod)
  
  
  # which cells in the functional responce matrix F[species,consumers] are non zero
  ind = which(results[[w]][[2]]$F[,] != 0, arr.ind = T)
  # create an empty (0) matrix with dimensions as the F matrix
  outflux = results[[w]][[2]]$F
  outflux[,] = 0
  
  if (nrow(ind)>0) { # this can occur if all plants are extinct
    for (i in 1:nrow(ind)) { # for all non-zero cells
      
      # calculate the energy that flows out of the resource (row) to the consumer (column)
      flux = results[[w]][[2]]$X[properties[[w]]$plant.rich + ind[i,][2]]*     # mass specific metabolic rate of consumer j
             results[[w]][[2]]$max_feed[ind[i,][2]]*                           # maximum feeding rate rel. to met. rate
             results[[w]][[1]][t.end,1+properties[[w]]$plant.rich+ind[i,][2]]* # biomass of cons. j at last time-step
             results[[w]][[2]]$F[ind[i,][1],ind[i,][2]]                        # functional responce at last time-step
      
      outflux[ind[i,][1],ind[i,][2]] = flux
      
    }
  } 
  colnames(outflux) = colnames(results[[w]][[1]])[-(1:(1+properties[[w]]$plant.rich))]
  rownames(outflux) = colnames(results[[w]][[1]])[-1]
  
  # calculate the energy that flows into each consumer (column) from each resource
  influx = results[[w]][[2]]$F
  influx[,] = 0
  
  if (nrow(ind)>0) { 
    for (i in 1:nrow(ind)) { 
      
      flux = outflux[ind[i,][1],ind[i,][2]]*
             results[[w]][[2]]$e[ind[i,][1]]  # prey i assimilation efficiency
      
      influx[ind[i,][1],ind[i,][2]] = flux
      
    }
  }
  colnames(influx) = colnames(results[[w]][[1]])[-(1:(1+properties[[w]]$plant.rich))]
  rownames(influx) = colnames(results[[w]][[1]])[-1]
  
  #,na.rm = T
  
  # herbivory control, measured as the ratio of outflows to inflows for herbivores
  properties[[w]]$herbivory.control = sum(outflux[grep("herb", rownames(outflux)), ])/
                                      sum( influx[ , grep("herb", colnames(influx))])
  
  # herbivory pressure, measured as total outflux from plants (incl. omnivores), per unit biomass
  properties[[w]]$herbivory.pressure = sum(outflux[1:properties[[w]]$plant.rich,])/ 
                                       properties[[w]]$plant.biomass
  
  # herbivore pressure, measured as total outflux from plants to herbivores, per unit biomass
  properties[[w]]$herbivore.pressure = sum(outflux[ , grep("herb", colnames(outflux))])/ 
                                       properties[[w]]$plant.biomass
  
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(w, '/', length(results)))
  #Sys.sleep(.05)
  if (w == length(results)) cat('- Done!')
}
beepr::beep(9)
properties = do.call(rbind.data.frame, properties)
properties[  seq(1, nrow(properties), 4),1] = "high for animals - high for plants"
properties[1+seq(1, nrow(properties), 4),1] = "high for animals - low for plants"
properties[2+seq(1, nrow(properties), 4),1] = "low for animals - high for plants"
properties[3+seq(1, nrow(properties), 4),1] = "low for animals - low for plants"
properties$scenario = as.factor(properties$scenario)
properties$ID = as.factor(rep(1:4000, each = 4))

saveRDS(properties, file="properties_20220721_N_40_comp_40000cont.RData")

properties1 = readRDS("properties_20220721_N_40_comp_40000.RData")
properties2 = readRDS("properties_20220721_N_40_comp_40000cont.RData")

properties1 = properties1[1:21618]
properties1 = do.call(rbind.data.frame, properties1)

properties2 = properties2[21619:24000]
properties2 = do.call(rbind.data.frame, properties2)

properties = rbind(properties1,properties2)

properties[  seq(1, nrow(properties), 6),1] = "high for animals - high inter for plants"
properties[1+seq(1, nrow(properties), 6),1] = "high for animals - high intra for plants"
properties[2+seq(1, nrow(properties), 6),1] = "high for animals - low for plants"
properties[3+seq(1, nrow(properties), 6),1] = "low for animals - high inter for plants"
properties[4+seq(1, nrow(properties), 6),1] = "low for animals - high intra for plants"
properties[5+seq(1, nrow(properties), 6),1] = "low for animals - low for plants"

properties$scenario = as.factor(properties$scenario)
properties$ID = as.factor(rep(1:4000, each = 6))

which(properties$plant.end==0)
length(which(properties$plant.end==0))
library(tidybayes)
library(see)
library(scico)
ggplot(properties, aes(x = log2(plant.rich), y=animal.biomass, color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=plant.end/plant.rich, 
           color = scenario)) +
  geom_point(aes(colour = scenario), 
             position = position_dodge(width = .9)) +
  stat_summary(fun = "mean",
               position = position_dodge(width = .9))
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()
  #geom_point(position = position_dodge(width = .9))

  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y=animals.end/40, 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9))  
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y=plant.end/plant.rich, 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern()
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y = log(omn.end/pred.end), 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern()
 
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y = log(omn.biomass/pred.biomass), 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern()
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = as.factor(log2(plant.rich)), 
             y = log(omn.end/pred.end), 
             color = scenario)) +
    geom_violinhalf(aes(color = scenario),
                    position = position_dodge(width = .8)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 1,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .8)) +
    theme_modern()
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y = log(herb.end/(omn.end+pred.end)), 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern() 
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y = log(herb.biomass/(omn.biomass+pred.biomass)), 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern()
  
  ggplot(properties[which(properties$plant.end!=0),], 
         aes(x = log2(plant.rich), 
             y = log10(herbivory.control), 
             color = scenario)) +
    geom_point(aes(fill = scenario), alpha = .1,
               position = position_jitterdodge(dodge.width = .9)) +
    geom_vline(xintercept = c(1.5,2.5,3.5), 
               linetype = "dashed") +
    stat_summary(fun = "mean",
                 geom = "point",
                 size = 4,
                 shape = 21,
                 aes(fill = scenario),
                 color = "black",
                 position = position_dodge(width = .9)) +
    theme_modern()
  
  
  
  
ggplot(properties, aes(x = log2(plant.rich), y=animals.end/60, color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.end), y=plant.biomass, color = scenario)) +
  geom_point(alpha = .1) +
  geom_smooth(method = "lm", formula = y ~ x) +
  #facet_wrap(~ scenario) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.end), y=primary.productivity, color = scenario)) +
  geom_point(alpha = .1) +
  geom_smooth(method = "lm", formula = y ~ x) +
  #facet_wrap(~ scenario) +
  theme_modern()


ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), y=animal.biomass, color = scenario)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2) ) +
  #facet_wrap(~ scenario) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=primary.productivity, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()


############################# PRIMARY PRODUCTIVITY #############################

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.end), 
           y = primary.productivity, 
           color = scenario)) +
  geom_point(alpha = .2, size = 3) +
  geom_smooth(aes(linetype = scenario),
              method = "lm", formula = y ~ x + I(x^2) ) +
  theme_modern() + 
  #facet_wrap(~ scenario) +
  scale_colour_manual(
    values = c("#893101","#893101","#EC9706","#EC9706")
  ) + 
  scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid"))

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.end), 
           y=primary.productivity, 
           color = scenario)) +
  geom_point(alpha = .2, size = 3) +
  geom_smooth(aes(linetype = scenario),
              method = "lm", formula = y ~ x + I(x^2)) +
  theme_modern() + 
  facet_wrap(~ scenario) +
  scale_colour_manual(
    values = c("#893101","#893101","#EC9706","#EC9706")
  ) + 
  scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid"))

################################################################################




ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log(herbivory.control + 1e-9), 
           y = (primary.productivity), 
           color = as.factor(plant.rich))) +
  geom_point(alpha = .5) +
  #geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~scenario) +
  #stat_pointinterval(aes(colour = scenario),
  #                   point_interval = "mean_qi",
  #                   position = position_dodge(width = .9)) +
  theme_modern()




# ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),], 
#        aes(x = log2(plant.rich), 
#            y=log10(herbivory.control + 1e-12), 
#            color = scenario)) +
#   #geom_point(position = position_dodge(width = .9)) +
#    stat_pointinterval(aes(colour = scenario),
#                       point_interval = "mean_qi",
#                       position = position_dodge(width = .9)) +
#   theme_modern()
# 
# ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),], 
#        aes(x = log2(plant.end), 
#            y=log10(herbivory.control + 1e-12), 
#            color = scenario)) +
#   geom_point(alpha = .5) +
#   geom_smooth(method = "lm") +
#   #stat_pointinterval(aes(colour = scenario),
#   #                   point_interval = "mean_qi",
#   #                   position = position_dodge(width = .9)) +
#   theme_modern()

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(pred.biomass/herb.biomass), 
           y = log10(herbivory.control + 1e-9), 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(y = "Herbivore control \n (Log-ratio of outflows/inflows for herbivores)",
       x = "Log-ratio of predator/herbivore biomass",
       color = "n plants")

##################### Predator-Herbivore ratio species #########################

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = plant.end, 
           y = log10(pred.end/herb.end), 
           color = plant.end)) +
  #stat_pointinterval(point_interval = "mean_qi") +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(x = "Number of plant species",
       y = "Log-ratio of predator/herbivore species",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0),] %>% arrange(plant.end), 
       aes(x = plant.end, 
           y = mean.predofherb, 
           color = plant.end)) +
  #stat_pointinterval(point_interval = "mean_qi") +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1,
                    begin = .05) +
  labs(x = "Number of plant species",
       y = "Average number of predators per herbivore",
       color = "n plants")


##################### Predator-Herbivore ratio biomass #########################

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = plant.end, 
           y = log10(pred.biomass/herb.biomass), 
           color = plant.end)) +
  #stat_pointinterval(point_interval = "mean_qi") +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(x = "Number of plant species",
       y = "Log-ratio of predator/herbivore biomass",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(pred.end/herb.end), 
           y = log10(pred.biomass/herb.biomass), 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(x = "Log-ratio of predator/herbivore species",
       y = "Log-ratio of predator/herbivore biomass",
       color = "n plants")

############### Net herbivore control ~ ratio of predators/herbivores ##########

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(pred.end/herb.end), 
           y = log10(herbivory.control + 1e-9), 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(y = "Herbivore control \n (Log-ratio of outflows/inflows for herbivores)",
       x = "Log-ratio of predator/herbivore species",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = mean.predofherb, 
           y = log10(herbivory.control + 1e-9), 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(#method = "lm", formula = y ~ x, 
              color = "#EC9706"
              ) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1, begin = .05) +
  labs(y = "Herbivore control \n (Log-ratio of outflows/inflows for herbivores)",
       x = "Average number of predators per herbivore",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log2(plant.rich), 
           y = log10(herbivory.control + 1e-9), 
           color = scenario)) +
  geom_point(alpha = .1, size = 3, position = position_jitter(.3)) +
  geom_smooth(method = "lm", formula = y ~ x #+ I(x^2)
    ,se = F#,color = "#EC9706"
  ) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~scenario) +
  theme_modern() + 
  #scale_color_scico(palette = 'lajolla',
  #                  direction = -1) +
  labs(y = "Herbivore control \n (Log-ratio of outflows/inflows for herbivores)",
       x = "Plant richness",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = (plant.end), 
           y=log10(herbivory.control + 1e-9), 
           color = scenario)) +
  geom_point(aes(fill = scenario), alpha = .1) +
  stat_summary(fun = "mean",
               geom = "point",
               size = 4,
               shape = 21,
               aes(fill = scenario),
               color = "black", alpha = .9,
               position = position_dodge(width = .1)) 


################### Herbivore pressure ~ herbivore biomass #####################

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(herb.biomass), 
           y = herbivore.pressure, 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  #geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(y = "Herbivore pressure \n (energy outflow  from plants to herbivores per unit plant biomass)",
       x = "Log herbivore biomass",
       color = "n plants")

################# Herbivore pressure ~ number of plant species #################

ggplot(properties[which(properties$plant.end!=0),] %>% arrange(plant.end), 
       aes(x = plant.end, 
           y = herbivore.pressure)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "#EC9706") +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~scenario) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1) +
  labs(y = "Herbivore pressure \n (energy outflow  from plants to herbivores per unit plant biomass)",
       x = "Number of plant species",
       color = "n plants")

ggplot(properties[which(properties$plant.end!=0),] %>% arrange(plant.end), 
       aes(x = plant.end, 
           y = herbivore.pressure, color = scenario)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~scenario) +
  theme_modern() + 
  #scale_color_scico(palette = 'lajolla',
  #                  direction = -1) +
  labs(y = "Herbivore pressure \n (energy outflow  from plants to herbivores per unit plant biomass)",
       x = "Number of plant species",
       color = "n plants")



ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=herbivory.pressure, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=herbivore.pressure, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = log2(plant.end), 
           y=herbivore.pressure, 
           color = scenario)) +
  geom_point(position = "jitter") +
  geom_smooth(method = "lm", 
              formula = y ~ x + I(x^2)) +
  
  #stat_pointinterval(aes(colour = scenario),
  #                   point_interval = "mean_qi",
  #                   position = position_dodge(width = .9)) +
  theme_modern()

ggplot(properties[which(properties$plant.end!=0 & properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(herbivory.control + 1e-9), 
           y=herbivore.pressure, 
           color = plant.end)) +
  geom_point(position = "jitter") +
  geom_smooth(method = "lm", 
              formula = y ~ x + I(x^2)) +
  facet_wrap(~scenario) +
  #stat_pointinterval(aes(colour = scenario),
  #                   point_interval = "mean_qi",
  #                   position = position_dodge(width = .9)) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1)

ggplot(properties[which(properties$plant.end!=0 & 
                          properties$herbivory.control!=Inf),] %>% arrange(plant.end), 
       aes(x = log10(herbivory.control + 1e-9), 
           y = herbivore.pressure, 
           color = plant.end)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", 
              formula = y ~ x + I(x^2), color = "#EC9706"
              ) +
  facet_wrap(~scenario) +
  #stat_pointinterval(aes(colour = scenario),
  #                   point_interval = "mean_qi",
  #                   position = position_dodge(width = .9)) +
  theme_modern() + 
  scale_color_scico(palette = 'lajolla',
                    direction = -1)

ggplot(properties[which(properties$plant.end!=0),], 
       aes(x = plant.biomass, 
           y = animal.biomass,
           color = scenario)) +
  geom_point(alpha = .7, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~scenario) +
  theme_modern() + 
  #scale_color_scico(palette = 'lajolla',
  #                  direction = -1) +
  labs(y = "animal.biomass",
       x = "plant.biomass",
       color = "n plants")

library(glmmTMB)
summary( 
  glmmTMB(herbivory.pressure ~
     log2(plant.rich) + 
     scenario +
     log2(plant.rich):scenario +
       (1|ID),
   data = properties[which(properties$plant.end!=0),])
)

summary( 
  glmmTMB(herbivory.control ~
            log2(plant.rich) + 
            scenario +
            log2(plant.rich):scenario +
            (1|ID),
          data = properties[which(properties$plant.end!=0),])
)


#summary( 
 m = glmmTMB(primary.productivity ~
            log2(plant.rich) + 
            scenario +
            log2(plant.rich):scenario +
            (1|ID),
          data = properties[which(properties$plant.end!=0),])
#)
library(performance)
check_model(m)





jac = Joacobian(results[[1]][[1]][601,-1],
                results[[1]][[2]]$ODE)



set.seed(123)
# first run a model to reach equilibrium
masses <- runif(20, 10, 100) #body mass of species
L <- create_Lmatrix(masses, 10, Ropt = 10)
L[L > 0] <- 1
mod <- create_model_Unscaled_nuts(20, 10, 3, masses, L)
mod <- initialise_default_Unscaled_nuts(mod, L)
biomasses <- masses ^ -0.75 * 10 ^ 4 #biomasses of species
biomasses <- append(runif(3, 20, 30), biomasses)
times <- seq(0, 100, 1)
sol <- lsoda_wrapper(times, biomasses, mod)
# get the final biomasses
final.bioms = sol[nrow(sol), -1]
# estimate jacobian
Joacobian(final.bioms, mod$ODE)