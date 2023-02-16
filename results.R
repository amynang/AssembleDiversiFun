library(tidyverse)


results = readRDS("results_2-16.RData")

properties = vector(mode = "list", length = length(results))

t.end = 201

for (w in 1:length(results)) {
  
  # turn biomass of extinct species to 0
  results[[w]][[1]][t.end,] = ifelse(results[[w]][[1]][t.end,]<=1e-6, 
                                     0, 
                                     results[[w]][[1]][t.end,])
  properties[[w]]$scenario = NA
  
  plant.rich = length(grep("plant", colnames(results[[w]][[1]])))
  properties[[w]]$plant.rich = plant.rich
  
  properties[[w]]$herb.start = length(grep("herb", colnames(results[[w]][[1]])))
  properties[[w]]$omn.start  = length(grep("omn",  colnames(results[[w]][[1]])))
  properties[[w]]$pred.start = length(grep("pred", colnames(results[[w]][[1]])))
  
  properties[[w]]$plant.end = length(which(results[[w]][[1]][t.end, grep("plant", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$animals.end = length(which(results[[w]][[1]][t.end, grep("herb|omn|pred", colnames(results[[w]][[1]]))]>0))
  
  properties[[w]]$herb.end = length(which(results[[w]][[1]][t.end, grep("herb", colnames(results[[w]][[1]]))]>0))
  properties[[w]]$omn.end  = length(which(results[[w]][[1]][t.end, grep("omn",  colnames(results[[w]][[1]]))]>0))
  properties[[w]]$pred.end = length(which(results[[w]][[1]][t.end, grep("pred", colnames(results[[w]][[1]]))]>0))
  
  
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
                  - results[[w]][[2]]$X[i]*
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

properties[  seq(1, nrow(properties), 6),1] = "high for animals - high inter for plants"
properties[1+seq(1, nrow(properties), 6),1] = "high for animals - high intra for plants"
properties[2+seq(1, nrow(properties), 6),1] = "high for animals - low for plants"
properties[3+seq(1, nrow(properties), 6),1] = "low for animals - high inter for plants"
properties[4+seq(1, nrow(properties), 6),1] = "low for animals - high intra for plants"
properties[5+seq(1, nrow(properties), 6),1] = "low for animals - low for plants"

properties$scenario = as.factor(properties$scenario)
properties$ID = as.factor(rep(1:4500, each = 6))

saveRDS(properties, file="properties_2-16.RData")

