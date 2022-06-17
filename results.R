library(tidyverse)


results = readRDS("results_20220609.RData")

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
for (w in 1:length(results)) {
  
  # turn biomass of extinct species to 0
  results[[w]][[1]][201,] = ifelse(results[[w]][[1]][201,]<=1e-6, 
                                   0, 
                                   results[[w]][[1]][201,])
  properties[[w]]$scenario = NA
  
  properties[[w]]$plant.rich = length(grep("plant", colnames(results[[w]][[1]])))
  
  properties[[w]]$plant.end = length(which(results[[w]][[1]][201, grep("plant", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$animals.end = length(which(results[[w]][[1]][201, grep("herb|omn|pred", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$herb.end = length(which(results[[w]][[1]][201, grep("herb", colnames(results[[w]][[1]]))]>0))
  properties[[w]]$omn.end  = length(which(results[[w]][[1]][201, grep("omn",  colnames(results[[w]][[1]]))]>0))
  properties[[w]]$pred.end = length(which(results[[w]][[1]][201, grep("pred", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$plant.biomass  = sum(results[[w]][[1]][201, grep("plant",         colnames(results[[w]][[1]]))])
  properties[[w]]$animal.biomass = sum(results[[w]][[1]][201, grep("herb|omn|pred", colnames(results[[w]][[1]]))])
  properties[[w]]$herb.biomass   = sum(results[[w]][[1]][201, grep("herb",          colnames(results[[w]][[1]]))])
  properties[[w]]$omn.biomass    = sum(results[[w]][[1]][201, grep("omn",           colnames(results[[w]][[1]]))])
  properties[[w]]$pred.biomass   = sum(results[[w]][[1]][201, grep("pred",          colnames(results[[w]][[1]]))])
  
  # calculate growth rate for each plant species
  plant.growth = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
  for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
    for (j in 1:dim(results[[w]][[2]]$alpha)[1]) {
      plant.growth[i] = 1 - (results[[w]][[2]]$alpha[i,j]* # competition from species j
                             results[[w]][[1]][201,1+j]/   # biomass of species j
                             results[[w]][[2]]$K)          # carrying capacity of species i
    }
  }
  #calculate plant productivity for each plant species
  plant.prod = vector(mode = "numeric", dim(results[[w]][[2]]$alpha)[1])
  for (i in 1:dim(results[[w]][[2]]$alpha)[1]) {
    plant.prod[i] = results[[w]][[2]]$r[i]*
                    plant.growth[i]*
                    results[[w]][[1]][201,1+i]
  }
  
  properties[[w]]$primary.productivity = sum(plant.prod)
  
  # calculate the energy that flows out of each resource (row) for each consumer (column)
  outflux = results[[w]][[2]]$F
  outflux[,] = NA
  for (i in 1:nrow(results[[w]][[2]]$F)) {
    for (j in 1:ncol(results[[w]][[2]]$F)) {
      
      outflux[i,j] = (results[[w]][[2]]$X[properties[[w]]$plant.rich + j]*   # mass specific metabolic rate of consumer j
                      results[[w]][[2]]$max_feed[j]*                         # maximum feeding rate rel. to met. rate
                      results[[w]][[1]][201,1+properties[[w]]$plant.rich+j]* # biomass of cons. j at last time-step
                      results[[w]][[2]]$F[i,j])/                             # functional responses at last timestep
                      results[[w]][[2]]$e[i]                                 # prey i assimilation efficiency
    }
  }
  colnames(outflux) = colnames(results[[w]][[1]])[-(1:(1+properties[[w]]$plant.rich))]
  rownames(outflux) = colnames(results[[w]][[1]])[-1]
  
  # calculate the energy that flows into each consumer (column) from each resource
  influx = results[[w]][[2]]$F
  influx[,] = NA
  for (i in 1:nrow(results[[w]][[2]]$F)) {
    for (j in 1:ncol(results[[w]][[2]]$F)) {
      
      influx[i,j] = (results[[w]][[2]]$X[properties[[w]]$plant.rich + j]*   # mass specific metabolic rate of consumer j
                     results[[w]][[2]]$max_feed[j]*                         # maximum feeding rate rel. to met. rate
                     results[[w]][[1]][201,1+properties[[w]]$plant.rich+j]* # biomass of cons. j at last time-step
                     results[[w]][[2]]$F[i,j])                              # functional responses at last timestep
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
  
}
beepr::beep(9)
porprties = do.call(rbind.data.frame, properties)
porprties[  seq(1, nrow(porprties), 4),1] = "early_early"
porprties[1+seq(1, nrow(porprties), 4),1] = "early_late"
porprties[2+seq(1, nrow(porprties), 4),1] = "late_early"
porprties[3+seq(1, nrow(porprties), 4),1] = "late_late"
porprties$scenario = as.factor(porprties$scenario)
porprties$ID = as.factor(rep(1:4000, each = 4))


which(porprties$plant.end==0)
length(which(porprties$plant.end==0))
library(tidybayes)
library(see)
ggplot(porprties, aes(x = log2(plant.rich), y=animal.biomass, color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(porprties, aes(x = log2(plant.rich), y=plant.end, color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     position = position_dodge(width = .9)) +
  theme_modern()
  #geom_point(position = position_dodge(width = .9))

ggplot(porprties[which(porprties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=primary.productivity, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()


ggplot(porprties[which(porprties$plant.end!=0 & porprties$herbivory.control!=Inf),], 
       aes(x = log2(plant.rich), 
           y=log10(herbivory.control + 1e-9), 
           color = scenario)) +
  #geom_point(position = position_dodge(width = .9)) +
   stat_pointinterval(aes(colour = scenario),
                      point_interval = "mean_qi",
                      position = position_dodge(width = .9)) +
  theme_modern()


ggplot(porprties[which(porprties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=herbivory.pressure, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()

ggplot(porprties[which(porprties$plant.end!=0),], 
       aes(x = log2(plant.rich), 
           y=herbivore.pressure, 
           color = scenario)) +
  stat_pointinterval(aes(colour = scenario),
                     point_interval = "mean_qi",
                     position = position_dodge(width = .9)) +
  theme_modern()

library(glmmTMB)
summary( 
  glmmTMB(herbivory.pressure ~
     log2(plant.rich) + 
     scenario +
     log2(plant.rich):scenario +
       (1|ID),
   data = porprties[which(porprties$plant.end!=0),])
)

summary( 
  glmmTMB(herbivory.control ~
            log2(plant.rich) + 
            scenario +
            log2(plant.rich):scenario +
            (1|ID),
          data = porprties[which(porprties$plant.end!=0),])
)


#summary( 
 m = glmmTMB(primary.productivity ~
            log2(plant.rich) + 
            scenario +
            log2(plant.rich):scenario +
            (1|ID),
          data = porprties[which(porprties$plant.end!=0),])
#)
library(performance)
check_model(m)
