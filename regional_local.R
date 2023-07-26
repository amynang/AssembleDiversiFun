library(tidyverse)
library(assembly)
library(igraph)
library(ATNr)
set.seed(321)

source("functions.R")

# the other plot foodweb function
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


# specify number of plants and plant consumers (herbivores & omnivores)
dim1 = 250
dim2 = 500

# generating nested plant-consumer interactions
# create a dataframe to generate nestedness
dum = data.frame(row = rep(1:dim1,dim2),
                 col = rep(1:dim2,each = dim1))
dum$strength = (dum$row/dim1)^2 + (dum$col/dim2)^2
# lattice::levelplot(strength ~ col*row, 
#                    data = dum)

# turn this into a matrix
duum = matrix(dum$strength,dim1,dim2)
# image(duum)

drum = duum
for (i in 1:(dim1*dim2)) {
  #probability of non-zero depends on the cell value, then value drawn from [.01,1]
  drum[i] = rbernoulli(1, (duum[i])^2)# * runif(1,.01,1)
}
# image(drum)

# shuffle them so that later on, this pattern will be unrelated to bodymasses
rand <- sample(ncol(drum))
#rand
drum = drum[,rand]

# create the interaction matrix
L = matrix(0,1000,1000)
# plants don't eat anything
L[,1:dim1] = 0

# randomly select 2/3 of animals to be plant consumers (incl. omnivores)
plant.cons = sample((dim1+1):1000, 500) %>% sort()

# put in the herbivory
L[1:dim1,plant.cons] = drum#[dim(drum)[1]:1, ]

n_species <- dim1+dim2+250
n_basal <- dim1

# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, -9, -3)),
                 sort(runif(n_species - n_basal, -9, 3)))


# create the allometric matrix
K <- create_Lmatrix(masses, 
                    n_basal, 
                    Ropt = 3.98,
                    gamma = 2, 
                    th = 0.01)

# predatory links are allometric
L[(dim1+1):n_species, ] = K[(dim1+1):n_species, ]

# remove cannibalism
diag(L) = 0


# split plant consumers to herbivores and omnivores
for (i in plant.cons) { 
  # the more generalist a plant consumer is, the more likely to be an omnivore
  L[(dim1+1):n_species, i] = L[(dim1+1):n_species, i]*(1-rbernoulli(1, sum(vegan::decostand(L[1:dim1,i],"pa"))/250))
}
# get the column indices of herbivores, omnivores, predators
herbiv = setdiff(which(colSums(L[(dim1+1):n_species, ]) == 0),1:dim1)
omniv  = setdiff(plant.cons, herbiv)
predat = setdiff((dim1+1):n_species, plant.cons)


#make the interaction matrix sparser (removing 30% of the interactions)
# the indices of non-zero cells
inds = which(L!=0)
# turn those cells to 0 with a probability 30%
L[inds] = as.integer(rbernoulli(length(inds),.7))*L[inds]
# check that consumers still have some resources
colSums(vegan::decostand(L[,(dim1+1):1000],"pa"))
min(colSums(vegan::decostand(L[,(dim1+1):1000],"pa")))
# marvel at the glory of your creation
#show_fw(vegan::decostand(L,"pa"))


# species names are 0000group000 where 0000 is the index of the species (relates to bodymass) group is
# plants,herbivores,omnivores,predators and 000 is the within group index
L = as.data.frame(L) %>%  rename_with(.cols = 1:dim1, ~c(paste0("_plant",sprintf("%03d", 1:250)))) %>% 
  rename_with(.cols = sort(herbiv), ~c(paste0("_herbiv",sprintf("%03d", 1:length(herbiv))))) %>% 
  rename_with(.cols = sort(predat), ~c(paste0("_predat",sprintf("%03d", 1:length(predat))))) %>% 
  rename_with(.cols = sort(omniv), ~c(paste0("_omniv",sprintf("%03d", 1:length(omniv))))) %>% 
  rename_with(.cols = 1:1000, ~c(paste0(sprintf("%04d", 1:1000), .)))
rownames(L) = colnames(L)

L = as.matrix(L)

#heatweb(L)



# create the 0/1 version of the food web
fw <- L
fw[fw > 0] <- 1


# an empty list
reg.loc = vector(mode = "list", 1)
# the first element of the list contains two objects:
# the interaction matrix of the regional meta-foodweb
# and the vector of bodymasses of the regional species
reg.loc[[1]][[1]] = L
reg.loc[[1]][[2]] = masses

# the subsequent elements will also contain two objects each:

# the first one is a local foodweb that is a random subset of the regional
# with 2-16 producers and 40 consumers
# consumers must have a resource and all producers have at least one consumer
# also, we want the subset to comprise a single foodweb (no isolated components)
# these will be our "early succession" foodwebs

# the second one will be the foodweb of a community with the same basal species
# and the same number of consumers as the first one but reduced similarity among
# consumers
# these will be our "late succession" foodwebs

# specify the levels of producer diversity
producer_diversity = seq(2,16,1)

# generate the "early succession" local foodwebs
# k iterations
for (k in 1:300) { 
  #for (j in consumer_diversity) {
  
  for (i in producer_diversity) {
    
    while(TRUE){ # while loop skips communities with isolated components
      # sample n=j local consumers from the regional pool
      local_consumers = sample(colnames(fw)[(n_basal+1):n_species], 40) %>% 
        sort()
      # sample n=i local producers from the regional pool
      local_producers = sample(colnames(fw)[1:n_basal], i) %>% 
        sort()
      # if the number of components is 1, we are good to proceed
      if(assembly:::.components(c(local_producers,local_consumers), fw) == 1 &
         length(assembly:::.find_isolated(c(local_producers,local_consumers), fw)) == 0
      ) break()
    }
    
    # subset the meta-fooweb to get the local foodweb
    local_fw = L[c(local_producers, local_consumers),
                 c(local_producers, local_consumers)]
    
    # add it to the list
    reg.loc[[length(reg.loc)+1]] = vector(mode = "list")
    reg.loc[[length(reg.loc)]][[1]] = local_fw
    
    #show_fw(vegan::decostand(local_fw,"pa"), title = "L-matrix model food web")
    
    #show_graph(colnames(local_fw), fw)
  }
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(k, '/', 4000))
  #Sys.sleep(.05)
  if (k == 4000) cat('- Done!')
}

# generate the "late succession" foodwebs
for (m in 2:length(reg.loc)) { 
  # each of the local foodwebs goes through a stochastic process that replaces 
  # species that are highly similar to other locals, with random species from the
  # regional pool, trying to reduce the overall linkage similarity
  sp_sim <- bride_of_similarity_filtering(colnames(reg.loc[[m]][[1]]), 
                                          fw, 
                                          t = .01, 
                                          max.iter = 500) %>% sort()
  
  # add the new foodweb next to the old one on the list
  reg.loc[[m]][[2]] = L[sp_sim, sp_sim]
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m-1, '/', length(reg.loc)-1))
  #Sys.sleep(.05)
  if (m == length(reg.loc)) cat('- Done!')
  
}


# for each early-late pair, add a matrix of high interspecific 
# a matrix of low interspecific but higher intraspecific competition
# and a matrix where interspecific competition decreases while intra- stays the same
set.seed(321)
for (m in 2:length(reg.loc)) {
  # number of plant species
  plants = length(grep("plant", colnames(reg.loc[[m]][[1]])))
  # add competition matrices
  # lower upper values specify **intraspecific** competition
  # interspecific for each plant then sums to 1-intraspesific
  reg.loc[[m]][[3]] = competition.N(lower = .5, upper = .6, plants) #high inter-
  reg.loc[[m]][[4]] = competition.N(lower = .9, upper = 1, plants) #low inter- high intra-
  reg.loc[[m]][[5]] = reg.loc[[m]][[4]]
  diag(reg.loc[[m]][[5]]) = diag(reg.loc[[m]][[3]]) # low overall
}


# save to working directory
saveRDS(reg.loc, file="reg.loc_2-16.RData")

