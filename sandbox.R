library(tidyverse)
library(assembly)
library(ATNr)
set.seed(321)
# 1) define number of species, their body masses, and the structure of the
# community
n_species <- 1000
n_basal <- 100
n_nut <- 16
# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
# 2) create the food web
# create the L matrix
L <- create_Lmatrix(masses, 
                    n_basal, 
                    Ropt = 3.98, #100, #3.98
                    gamma = 2, 
                    th = 0.01)
diag(L) = 0

N  <-  length(L)/2                              # the number of random values to replace
inds <- round ( runif(N, 1, length(L)) )   # draw random values from [1, length(L)]
L[inds] <- 0                               # use the random values as indicies to L, for which to replace





colnames(L) = c(paste0("broducer",sprintf("%04d", 1:n_basal)),
                paste0("consumer",sprintf("%04d", 1:(n_species-n_basal)))
                )

rownames(L) = colnames(L)

# create the 0/1 version of the food web
fw <- L
fw[fw > 0] <- 1

show_fw(fw, title = "L-matrix model food web")
# 
# S = 50
# sp <- draw_random_species(S, colnames(fw)[17:516])
# 
# sp = 
# sp_resource <- resource_filtering(sp, adirondack, keep.n.basal = TRUE)

# local_consumers = sample(colnames(fw)[(n_basal+1):n_species], 50) %>% sort()
# 
# local_fw = fw[c(colnames(fw)[1:n_basal], local_consumers),
#               c(colnames(fw)[1:n_basal], local_consumers)]
# 
# show_fw(local_fw, title = "L-matrix model food web")
# 
# 
# show_graph(colnames(local_fw), fw)




# sp_sim <- similarity_filtering(colnames(local_fw), 
#                                fw, 
#                                t = 0, 
#                                max.iter = 100)
# 
# 
# show_graph(sp_sim, fw)


# 4 levels of producer diversity
producer_diversity = c(2,4,8,16)
# 4 levels of consumer diversity
#consumer_diversity = c(30,40,50,60)
# an empty list 
local_fws = vector(mode = "list")

# k iterations
for (k in 1:1000) { 
  #for (j in consumer_diversity) {
    
    for (i in producer_diversity) {
      
      while(TRUE){ # while loop skips communities with isolated components
        # sample n=j local consumers from the regional pool
        local_consumers = sample(colnames(fw)[(n_basal+1):n_species], 60) %>% 
          sort()
        # sample n=i local producers from the regional pool
        local_producers = sample(colnames(fw)[1:n_basal], i) %>% 
          sort()
        # if the number of components is 1, we are good to proceed
        if(assembly:::.components(c(local_producers,local_consumers), fw) == 1) break()
      }
      

      # subset the meta-fooweb to get the local foodweb
      local_fw = L[c(local_producers, local_consumers),
                   c(local_producers, local_consumers)]
      #diag(local_fw) = 0
      # Do I need this?
      # try(
      #sp_resource <- resource_filtering(local_fw, fw, keep.n.basal = TRUE)
      #   )
      # 
      # local_fw = fw[sp_resource,
      #               sp_resource]
      
      # add that foodweb to the list
      local_fws[[length(local_fws) + 1]] = local_fw
      
      show_fw(vegan::decostand(local_fw,"pa"), title = "L-matrix model food web")
      
      #show_graph(colnames(local_fw), fw)
    }
  #}
}
gg = vector(mode = "numeric", length=length(local_fws))
for (g in 1:length(local_fws)) {
  gg[g] = assembly:::.components(colnames(local_fws[[g]]), fw)
}
length(which(gg>1))
# loc_fws = local_fws[-which(gg>1)]
# View(local_fws[-which(gg>1)])

late_succession = vector(mode = "list", length = length(local_fws))

# for (m in 1:length(late_succession)) { 
#   
#   sp_sim <- similarity_filtering(colnames(local_fws[[m]]), 
#                                           fw, 
#                                           t = 0, 
#                                           max.iter = 100) %>% 
#     sort()
#   
#   late_succession[[m]] <- L[sp_sim,
#                             sp_sim]
#   
#   
#   cat('\014')
#   #cat(paste0(round((m/1600)*100), '%'))
#   cat(paste0(m, '/160'))
#   #Sys.sleep(.05)
#   if (m == length(local_fws)) cat('- Done!')
#   
# }

for (m in 1:length(late_succession)) { 
  
  # this will give 0x0 matrices in case of error (isolated species or components)
  sp_sim <- tryCatch(similarity_filtering(colnames(local_fws[[m]]), 
                                          fw, 
                                          t = 0, 
                                          max.iter = 1000) %>% 
                       sort(),
                     error = function(e) 0)
  
  late_succession[[m]] <- L[sp_sim,
                            sp_sim]
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m, '/', length(late_succession)))
  #Sys.sleep(.05)
  if (m == length(local_fws)) cat('- Done!')
  
}

beepr::beep(9)


show_fw(vegan::decostand(local_fws[[15]],"pa"), title = "L-matrix model food web")
show_fw(vegan::decostand(late_succession[[15]],"pa"), title = "L-matrix model food web")


# saveRDS(late_succession, file="late_succession.RData")
# saveRDS(local_fws, file="local_fws.RData")
# 
# late_succession = readRDS("late_succession.RData")

dims = vector(mode = "numeric", length(late_succession))

for (i in 1:length(late_succession)) {
  dims[i] = dim(late_succession[[i]])[1]
}
  

model_unscaled_nuts <- create_model_Unscaled_nuts(62, 2, 2, 
                                                  masses[which(colnames(L) %in% colnames(local_fws[[1]]))], 
                                                  vegan::decostand(local_fws[[1]],"pa"))

# for a model created by create_model_Unscaled_nuts():
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, 
                                                        local_fws[[1]])
model_unscaled_nuts$initialisations()


biomasses <- masses[which(colnames(L) %in% colnames(local_fws[[1]]))] ^ -0.75 * 1e1 # starting biomasses
biomasses <- append(runif(2, 20, 30), biomasses) # nutrient concentration
# defining the desired integration time
times <- seq(0, 1500, 1)
sol <- lsoda_wrapper(times, biomasses, model_unscaled_nuts)
plot_odeweb(sol, model_unscaled_nuts$nb_s)

colnames(sol) = c("time",
                   paste0("nut",1:2),
                   paste0("plant",1:2),
                   paste0("animal",1:60))

dat = as.data.frame(sol) %>% pivot_longer(!time, 
                                          names_to = "species", 
                                          values_to = "biomass")
dat$taxon = substr(dat$species, 1, 3)


ggplot2::ggplot(dat, aes(x=time, 
                         y=biomass, 
                         color = taxon)) +
  geom_point()




model_unscaled_nuts <- create_model_Unscaled_nuts(62, 2, 2, 
                                                  masses[which(colnames(L) %in% colnames(late_succession[[1]]))], 
                                                  vegan::decostand(late_succession[[1]],"pa"))

# for a model created by create_model_Unscaled_nuts():
model_unscaled_nuts <- initialise_default_Unscaled_nuts(model_unscaled_nuts, 
                                                        late_succession[[1]])
model_unscaled_nuts$initialisations()


biomasses <- masses[which(colnames(L) %in% colnames(late_succession[[1]]))] ^ -0.75 * 1e1 # starting biomasses
biomasses <- append(runif(2, 20, 30), biomasses) # nutrient concentration
# defining the desired integration time
times <- seq(0, 1500, 1)
sol <- lsoda_wrapper(times, biomasses, model_unscaled_nuts)
plot_odeweb(sol, model_unscaled_nuts$nb_s)

colnames(sol) = c("time",
                  paste0("nut",1:2),
                  paste0("plant",1:2),
                  paste0("animal",1:60))

dat = as.data.frame(sol) %>% pivot_longer(!time, 
                                          names_to = "species", 
                                          values_to = "biomass")
dat$taxon = substr(dat$species, 1, 3)


ggplot2::ggplot(dat, aes(x=time, 
                         y=biomass, 
                         color = taxon)) +
  geom_point()




# running simulations for the Schneider model
sol <- deSolve::lsoda(
  biomasses,
  times,
  function(t, y, params) {
    return(list(params$ODE(y, t)))
  },
  model_unscaled_nuts
)

sp_sim <- similarity_filtering(colnames(local_fws[[65]]),
                               fw,
                               t = 0,
                               max.iter = 100) %>% 
  sort()


show_fw(local_fws[[65]], title = "L-matrix model food web")
show_fw(fw[sp_sim,
           sp_sim], title = "L-matrix model food web")


sp_sim <- bride_of_similarity_filtering(colnames(local_fws[[5]]),
                                        fw,
                                        t = 0,
                                        max.iter = 500) %>% 
  sort()


show_fw(local_fws[[5]], title = "L-matrix model food web")
show_fw(fw[sp_sim,
           sp_sim], title = "L-matrix model food web")




bride_of_similarity_filtering(colnames(local_fws[[1]]), 
                              fw, 
                              t = 0, 
                              max.iter = 1000)





# function to plot the fw
show_fw <- function(mat, title = NULL) {
  par(mar = c(.5, .5, 2, .5))
  S <- nrow(mat)
  mat <- mat[nrow(mat):1, ]
  mat <- t(mat)
  image(mat, col = c("goldenrod", "steelblue"),
        frame = FALSE, axes = FALSE)
  title(title)
  grid(nx = S, ny = S, lty = 1, col = adjustcolor("grey20", alpha.f = .1))
}


bride_of_similarity_filtering <- function (sp.names, metaweb, t = 0, method = "jaccard", max.iter = 1000) 
{
  new_sp <- sp.names
  
  while(TRUE){
    for (i in seq_len(max.iter)) new_sp <- assembly:::.move(new_sp, metaweb, 
                                                            t = t)
    if(assembly:::.components(new_sp, metaweb) == 1) break()
  }
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  # if (.components(new_sp, metaweb) > 1) 
  #   stop("Isolated component detected")
  return(new_sp)
}


bride_of_resource_filtering <- function (sp.names, metaweb, keep.n.basal = FALSE) 
{
  isolated <- assembly:::.find_isolated(sp.names, metaweb)
  new_sp <- sp.names
  i <- 0
  while (length(isolated) != 0 & i < 1e+06) {
    replacements <- assembly:::.find_replacements(new_sp, isolated, 
                                       metaweb, keep.n.basal)
    new_sp <- union(setdiff(new_sp, isolated), replacements)
    isolated <- assembly:::.find_isolated(new_sp, metaweb)
    i <- i + 1
  }
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  # if (assembly:::.components(new_sp, metaweb) != 1) 
  #   stop("Isolated component detected")
  return(new_sp)
}


dum = matrix(1, 
             ncol = 10, 
             nrow = 10)

colnames(dum) = c(paste0("broducer",sprintf("%04d", 1:4)),
                paste0("consumer",sprintf("%04d", 1:6))
)

rownames(dum) = colnames(dum)


local_consumers = sample(colnames(dum)[(4+1):10], 3) %>% 
  sort()
# sample n=i local producers from the regional pool
local_producers = sample(colnames(dum)[1:4], 2) %>% 
  sort()

local_fw = fw[c(local_producers, local_consumers),
              c(local_producers, local_consumers)]
sp_resource <- resource_filtering(local_fw, dum, keep.n.basal = TRUE)



##### Trait based assembly #####
sp = 5
tr = 3
trait.mat = matrix(rnorm(sp*tr,0,1), sp, tr, byrow = T)
hist(trait.mat[,1])

vegan::vegdist(trait.mat)
ecodist::distance(trait.mat)



# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
masses <- sort(rtgamma(1e3, 1, scale=2, a=0, b=10))^3

# 2) create the food web
# create the L matrix
L <- create_Lmatrix(masses, 
                    0, 
                    Ropt = 3.98, #100, #3.98
                    gamma = 2, 
                    th = 0)

L <- create_matrix(masses, 
                   0, 
                   Ropt = 3.98, #100, #3.98
                   gamma = 2, 
                   th = 0)

rownames(L) = colnames(L) = paste0("sp",1:length(masses))

#show_fw(vegan::decostand(L,"pa"))
heatweb(L)


heatweb <- function(mat) {
  heat <- mat %>% #na_if(., 0) %>% 
    as.data.frame() %>%
    rownames_to_column("id") %>%
    pivot_longer(-c(id), names_to = "species", values_to = "strength") %>%
    mutate(species= fct_relevel(species,colnames(mat))) %>%
    ggplot(aes(x=species, y=ordered(id, levels = rev(unique(id))), fill=strength)) + 
    geom_raster() +
    #theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none") +
    #scico::scale_fill_scico(palette = "lajolla")
    #scale_fill_distiller(palette = "Spectral", direction = -1)
    scale_fill_viridis_c(option = "inferno", direction = -1)
  return(heat)
}

create_matrix <- function (BM, nb_b, Ropt = 100, gamma = 2, th = 0.01) 
{
  Lmatrix <- function(BM, nb_b, Ropt, gamma, th) {
    s <- length(BM)
    L <- matrix(rep(BM, s), s, s, byrow = TRUE)/(matrix(rep(BM, 
                                                            s), s, s) * Ropt)
    L <- (L * exp(1 - L))^gamma
    L[L < th] <- 0
    L[, 1:nb_b] <- 0
    return(L)
  }
  s <- length(BM)
  L <- Lmatrix(BM, nb_b, Ropt, gamma, th)
  isolated <- ifelse(any(colSums(L) + rowSums(L) == 0), TRUE, 
                     FALSE)
  cons_no_prey <- ifelse(any(colSums(L[, (nb_b + 1):s]) == 
                               0), TRUE, FALSE)
  tro_lev <- tryCatch(ATNr::TroLev(L), error = function(e) NULL)
  connected <- is_connected(graph_from_adjacency_matrix(vegan::decostand(L,"pa")))
  # i <- 0
  # while ((isolated | cons_no_prey | is.null(tro_lev) | !connected) & 
  #        i < 100) {
  #   L <- Lmatrix(BM, nb_b, Ropt, gamma, th)
  #   isolated <- ifelse(any(colSums(L) + rowSums(L) == 0), 
  #                      TRUE, FALSE)
  #   cons_no_prey <- ifelse(any(colSums(L[, (nb_b + 1):s]) == 
  #                                0), TRUE, FALSE)
  #   if (!isolated) {
  #     tro_lev <- tryCatch(ATNr::TroLev(L), error = function(e) NULL)
  #   }
  #   i <- i + 1
  # }
  if (isolated) 
    warning("Presence of an isolated species after 100 iterations.")
  if (cons_no_prey) 
    warning("Presence of consumer without prey after 100 iterations.")
  if (!connected) 
    warning("several conected component detected")
  return(L)
}



n = 3

set.seed(321)
biomasses <- runif(dim(reg.loc[[n]][[2]])[1], 2, 3) # starting biomasses

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
model_scaled$ext = 0.001
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
# solll[solll==""] = "nut"
ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
  geom_point()




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
# solll[solll==""] = "nut"
ggplot2::ggplot(solll[,], aes(x=time, y=biomass, color = taxon)) +
  geom_point()







n=2

model_unscaled_nuts <- create_model_Unscaled_nuts(# number of species
  dim(reg.loc[[n]][[2]])[1], 
  # number of basal species
  length(grep("plant", colnames(reg.loc[[n]][[2]]))), 
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

