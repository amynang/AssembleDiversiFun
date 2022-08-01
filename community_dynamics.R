library(tidyverse)
#library(assembly)
library(ATNr)
set.seed(321)

reg.loc = readRDS("reg.loc_20220721_N_40_comp.RData")

# the first element of the list contains two objects:
# the interaction matrix of the regional meta-foodweb
# and the vector of bodymasses of the regional species

# the subsequent elements contain four objects each:

# the first one is a local foodweb that is a random subset of the regional
# with 2, 4, 8 or 16 producers and 60 consumers
# consumers must have a resource and all producers have at least one consumer
# also, we want the subset to comprise a single foodweb (no isolated components)
# these are our "early succession" foodwebs

# the second one is the foodweb of a community with the same basal species
# and the same number of consumers as the first one but reduced similarity among
# consumers
# these are our "late succession" foodwebs

# the third and fourth ones are plant competition matrices
# with high interspecific low intraspecific and low interspecific high intraspecific
# competition respactively, such that the overlall competition that a species 
# experiences (the column sum) remains 1

results = vector(mode = "list")

#results[[1]][[1]] = 1

for (i in 2:length(reg.loc)) { # for each food-web
  for (j in 1:2) { # random or dissimilar consumers
    for (k in 3:5) { # high or low interspecific plant competition
      # we draw random starting biomasses from U(2,3)
      biomasses <- runif(dim(reg.loc[[i]][[j]])[1], .02, .03)
      # number of species
      species = dim(reg.loc[[i]][[j]])[1]
      # number of basal species
      plants = length(grep("plant", colnames(reg.loc[[i]][[j]])))
      # bodymasses of species
      bodymasses = reg.loc[[1]][[2]][as.integer(substr(colnames(reg.loc[[i]][[j]]), 1, 4))]
      
      # like Delmas
      model_scaled <- create_model_Scaled(# number of species
                                          species, 
                                          # number of basal species
                                          plants, 
                                          # bodymasses of species
                                          bodymasses, 
                                          # binary interaction matrix
                                          vegan::decostand(reg.loc[[i]][[j]],"pa")
      )
      
      # set.seed(i) # so that food-webs in every group differ only in animal/plant competition
      model_scaled <- initialise_default_Scaled(model_scaled)
      model_scaled$initialisations()
      # change W so that unlike Delmas, generalists are as efficient as specialists
      #model_scaled$w = vegan::decostand(model_scaled$w,"pa")
      
      # reduce interference competition
      #model_scaled$c = rep(.3,60)
      
      # plant competition scenario
      model_scaled$alpha = reg.loc[[i]][[k]]
      
      # timesteps
      times <- seq(0, 40000, by = 200)
      # solve
      sol <- lsoda_wrapper(times, biomasses, model_scaled, verbose = FALSE)
      soll = as.data.frame(sol)
      
      colnames(soll) = c("time", colnames(reg.loc[[i]][[j]]))
      
      # for each foodweb save biomass changes and parameters at last timestep
      results[[length(results)+1]] = vector(mode = "list")
      results[[length(results)]][[1]] = soll
      results[[length(results)]][[2]] = as.list(model_scaled)
      
      #print(sum(soll[201,-1]))
      #print(max(Re(eigen(Joacobian(sol[nrow(sol), -1], model_scaled$ODE))$values)))
    }
  }
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(i, '/', length(reg.loc)-1))
  #Sys.sleep(.05)
  if (i == length(reg.loc)-1) cat('- Done!')  # if (i == length(reg.loc)-1) cat('- Done!')
}


# save to working directory
saveRDS(results, file="results_20220721_N_40_comp_40000.RData")

# for (i in 2:ncol(results[[1]][[1]])) {
#   
#   
# }
# 
# colnames(soll) = c("time",
#                    # paste0("nut_",1:2),
#                    colnames(reg.loc[[n]][[s]]))
# 
# solll = soll %>% pivot_longer(!time, names_to = "species", values_to = "biomass")
# solll$taxon = substr(solll$species, 6, 8)
# 
# # plot results
# ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
#   geom_point() 
