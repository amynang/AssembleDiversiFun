library(tidyverse)
library(assembly)
library(ATNr)
set.seed(321)

reg.loc = readRDS("reg.loc_20220506.RData")
# the first element of the list contains two objects:
# the interaction matrix of the regional meta-foodweb
# and the vector of bodymasses of the regional species

# the subsequent elements also contain two objects each:

# the first one is a local foodweb that is a random subset of the regional
# with 2, 4, 8 or 16 producers and 60 consumers
# consumers must have a resource and all producers have at least one consumer
# also, we want the subset to comprise a single foodweb (no isolated components)
# these are our "early succession" foodwebs

# the second one is the foodweb of a community with the same basal species
# and the same number of consumers as the first one but reduced similarity among
# consumers
# these are our "late succession" foodwebs



# Pick a number, any number (1,41]
n = 3

set.seed(321)
biomasses <- runif(dim(reg.loc[[n]][[2]])[1], 2, 3) # starting biomasses

# all examples below are with reg.loc[[n]][[2]] so the late succession version 
# of the chosen community

# like Delmas (but w/ allometric rel. preference matrix)
model_scaled <- create_model_Scaled(# number of species
                                    dim(reg.loc[[n]][[2]])[1], 
                                    # number of basal species
                                    length(grep("plant", colnames(reg.loc[[n]][[2]]))), 
                                    # bodymasses of species
                                    reg.loc[[1]][[2]][as.integer(substr(colnames(reg.loc[[n]][[2]]), 1, 4))], 
                                    # binary interaction matrix
                                    vegan::decostand(reg.loc[[n]][[2]],"pa")
)

model_scaled <- initialise_default_Scaled(model_scaled)
model_scaled$initialisations()

# change W so that unlike Delmas, relative preferences are allometric
model_scaled$w = vegan::decostand(reg.loc[[n]][[2]][,-grep("plant", colnames(reg.loc[[n]][[2]]))], "total" , 2)
str(model_scaled)

#model_scaled$alpha = matrix(.4,4,4)
#diag(model_scaled$alpha) = 1

times <- seq(0, 1000, by = 5)
#biomasses <- runif(dim(reg.loc[[n]][[2]])[1], 2, 3) # starting biomasses

sol1 <- lsoda_wrapper(times, biomasses, model_scaled, verbose = FALSE)
#plot_odeweb(sol1, dim(local_fws[[n]])[1])
soll = as.data.frame(sol1)
colnames(soll) = c("time",
                   # paste0("nut_",1:2),
                   colnames(reg.loc[[n]][[2]]))

solll = soll %>% pivot_longer(!time, names_to = "species", values_to = "biomass")
solll$taxon = substr(solll$species, 6, 8)

# plot results
ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
  geom_point()



# like Binzer
model_unscaled <- create_model_Unscaled(# number of species
                                        dim(reg.loc[[n]][[2]])[1], 
                                        # number of basal species
                                        length(grep("plant", colnames(reg.loc[[n]][[2]]))), 
                                        # bodymasses of species
                                        reg.loc[[1]][[2]][as.integer(substr(colnames(reg.loc[[n]][[2]]), 1, 4))], 
                                        # binary interaction matrix
                                        vegan::decostand(reg.loc[[n]][[2]],"pa")
                                        )

model_unscaled <- initialise_default_Unscaled(model_unscaled)
model_unscaled$initialisations()
str(model_unscaled)

times <- seq(0, 1000, by = 5)
#biomasses <- runif(dim(reg.loc[[n]][[2]])[1], 2, 3) # starting biomasses

sol1 <- lsoda_wrapper(times, biomasses, model_unscaled, verbose = FALSE)
#plot_odeweb(sol1, dim(local_fws[[n]])[1])
soll = as.data.frame(sol1)
colnames(soll) = c("time",
                   # paste0("nut_",1:2),
                   colnames(reg.loc[[n]][[2]]))

solll = soll %>% pivot_longer(!time, names_to = "species", values_to = "biomass")
solll$taxon = substr(solll$species, 6, 8)

# plot results
ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
  geom_point()


# like Schneider 
model_unscaled_nuts <- create_model_Unscaled_nuts(# number of species
  dim(reg.loc[[n]][[2]])[1], 
  # number of basal species
  length(grep("plant", colnames(reg.loc[[n]][[2]]))), 
  # number of nutrients
  4,
  # bodymasses of species
  reg.loc[[1]][[2]][as.integer(substr(colnames(reg.loc[[n]][[2]]), 1, 4))], 
  # binary interaction matrix
  vegan::decostand(reg.loc[[n]][[2]],"pa")
)

model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, reg.loc[[n]][[2]])
model_unscaled_nuts$initialisations()
str(model_unscaled_nuts)

times <- seq(0, 1000, by = 5)
set.seed(321)
biomasses <- runif(dim(reg.loc[[n]][[2]])[1], 2, 3) # starting biomasses
biomasses <- append(runif(4, 2, 3), biomasses) # nutrient concentration

sol1 <- lsoda_wrapper(times, biomasses, model_unscaled_nuts, verbose = FALSE)
#plot_odeweb(sol1, dim(local_fws[[n]])[1])
soll = as.data.frame(sol1)
colnames(soll) = c("time",
                   paste0("nut_",1:4),
                   colnames(reg.loc[[n]][[2]]))

solll = soll %>% pivot_longer(!time, names_to = "species", values_to = "biomass")
solll$taxon = substr(solll$species, 6, 8)
solll[solll==""] = "nut"
ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
  geom_point()

