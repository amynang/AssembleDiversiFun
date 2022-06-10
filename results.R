
w = 16
outflux = results[[w]][[2]]$F#[i,j]
outflux[,] = NA
for (i in 1:nrow(results[[w]][[2]]$F)) {
  for (j in 1:ncol(results[[w]][[2]]$F)) {
    
    outflux[i,j] = (results[[w]][[2]]$X[j]* #mass specific metabolic rate
                    results[[w]][[2]]$max_feed[16+j]* # maximum feeding rate rel. to met. rate
                    results[[w]][[1]][201,16+j]* # biomass at last time-step
                    results[[w]][[2]]$F[i,j])/ # functional responses at last timestep
                    results[[w]][[2]]$e[j] # assimilation efficiency
  }
}
colnames(outflux) = colnames(results[[w]][[1]])[-(1:17)]
rownames(outflux) = colnames(results[[w]][[1]])[-1]
which(colSums(outflux)==0)
sum(outflux[1:16,],na.rm = T)


results[[w]][[1]][201,2:3]
View(results[[w]][[1]])

i=16
properties = vector(mode = "list", length = length(results))
for (i in 1:length(results)) {
  properties[[i]]$scenario = NA
  
  properties[[i]]$plant.rich = length(grep("plant", colnames(results[[i]][[1]])))
  
  properties[[i]]$plant.end = length(which(results[[i]][[1]][201, grep("plant", colnames(results[[i]][[1]]))]>1e-6))
  
  properties[[i]]$animals.end = length(which(results[[i]][[1]][201, grep("herb|omn|pred", colnames(results[[i]][[1]]))]>1e-6))
  
  properties[[i]]$plant.biomass = sum(results[[i]][[1]][201, 1+which(results[[i]][[1]][201, grep("plant", colnames(results[[i]][[1]]))]>1e-6)])
  
  properties[[i]]$animal.biomass = sum(results[[i]][[1]][201, 1 + length(grep("plant", colnames(results[[i]][[1]]))) +
                                                           which(results[[i]][[1]][201, grep("herb|omn|pred", colnames(results[[i]][[1]]))]>1e-6)])
  
}
porprties = do.call(rbind.data.frame, properties)
porprties[  seq(1, nrow(porprties), 4),1] = "early_early"
porprties[1+seq(1, nrow(porprties), 4),1] = "early_late"
porprties[2+seq(1, nrow(porprties), 4),1] = "late_early"
porprties[3+seq(1, nrow(porprties), 4),1] = "late_late"
porprties$scenario = as.factor(porprties$scenario)


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
