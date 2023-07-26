library(tidyverse)
#library(assembly)
library(ATNr)
set.seed(321)



# unlike the default option in the scaled model, we want to scale all rates rel.
# to those of the smallest producer *in the regional pool*
# this way we can rate-based attributes across food-webs
initialise_default_Scaled.local <- function (model, min.BM) 
{
  utils::data("schneider", envir = environment())
  schneider[["nb_s"]] <- model$nb_s
  schneider[["nb_b"]] <- model$nb_b
  schneider[["BM"]] <- model$BM
  ar <- 1
  K <- 10
  w <- sweep(x = model$fw, MARGIN = 2, FUN = "/", colSums(model$fw))
  model$w <- w[, (model$nb_b + 1):model$nb_s]
  #min.BM = min(model$BM[1:model$nb_b])
  # min.BM <- with(schneider, min(BM[1:nb_b]))
  model$X <- with(schneider, 0.314 * BM^-0.25 / min.BM^-0.25)
  # model$X[1:schneider$nb_b] <- 0.0
  model$e <- with(schneider, c(rep(e_P, nb_b), rep(e_A, nb_s - 
                                                     nb_b)))
  model$c <- rep(0.8, model$nb_s - model$nb_b)
  model$max_feed <- rep(8, model$nb_s - model$nb_b)
  model$B0 <- rep(0.5, model$nb_s - model$nb_b)
  model$q <- stats::rnorm(1, 1.2, 0.2)
  model$r <- with(schneider, (ar * BM[1:nb_b]^-0.25)/(ar * 
                                                        min.BM^-0.25))
  model$K <- 10
  model$F <- with(schneider, matrix(0, nrow = model$nb_s, 
                                    ncol = model$nb_s - model$nb_b))
  model$alpha <- matrix(0, nrow = model$nb_b, ncol = model$nb_b)
  diag(model$alpha) = 1
  return(model)
}


reg.loc = readRDS("reg.loc_2-16.RData")
# the first element of the list contains two objects:
# the interaction matrix of the regional meta-foodweb
# and the vector of bodymasses of the regional species

# the subsequent elements contain five objects each:

# the first one is a local foodweb that is a random subset of the regional
# with 2 to 16 producers and 40 consumers
# consumers must have a resource and all producers have at least one consumer
# also, we want the subset to comprise a single foodweb (no isolated components)
# these are our "early succession" foodwebs

# the second one is the foodweb of a community with the same basal species
# and the same number of consumers as the first one but reduced similarity among
# consumers
# these are our "late succession" foodwebs

# the third, fourth and fifth ones are plant competition matrices
# with high interspecific low intraspecific, low interspecific high intraspecific
# and low inter and intra-
# competition respectively


results = vector(mode = "list")

for (i in 2:length(reg.loc)) { # for each food-web
  for (j in 1:2) { # random or dissimilar consumers
    for (k in 3:5) { # high or low interspecific plant competition
      
      min.BM = reg.loc[[1]][[2]][1]
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
      
      model_scaled <- initialise_default_Scaled.local(model_scaled, min.BM = min.BM)
      model_scaled$initialisations()
      
      # fix the Hill exponent
      model_scaled$q <- 1.2
      
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
saveRDS(results, file="results_2-16.RData")

