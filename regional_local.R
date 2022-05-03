library(tidyverse)
library(assembly)
library(ATNr)
set.seed(123)

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

# quick & dirty way for nested plant-consumer interactions
# create a dataframe to generate nestedness
dum = data.frame(row = rep(1:dim1,dim2),
                 col = rep(1:dim2,each = dim1))
dum$strength = (dum$row/dim1)^2 + (dum$col/dim2)^2
lattice::levelplot(strength ~ col*row, 
                   data = dum)

# turn this into a matrix
duum = matrix(dum$strength,dim1,dim2)
image(duum)

drum = duum
for (i in 1:(dim1*dim2)) {
  #probability of non-zero depends on the cell value, then value drawn from [.01,1]
  drum[i] = rbernoulli(1, (duum[i])^2) * runif(1,.01,1)
}
image(drum)

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
# masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
#                  sort(runif(n_species - n_basal, 2, 9)))

# from 1ngram to 1kgram
masses <- 10 ^ c(sort(runif(n_basal, -9, -3)),
                 sort(runif(n_species - n_basal, -9, 3)))


# create the allometric matrix
K <- create_Lmatrix(masses, 
                    n_basal, 
                    Ropt = 3.98, #100, #3.98
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
show_fw(vegan::decostand(L,"pa"))


# species names are 0000group000 where 0000 is the index of the species (relates to bodymass) group is
# plants,herbivores,omnivores,predators and 000 is the within group index
L = as.data.frame(L) %>%  rename_with(.cols = 1:dim1, ~c(paste0("_plant",sprintf("%03d", 1:250)))) %>% 
  rename_with(.cols = sort(herbiv), ~c(paste0("_herbiv",sprintf("%03d", 1:length(herbiv))))) %>% 
  rename_with(.cols = sort(predat), ~c(paste0("_predat",sprintf("%03d", 1:length(predat))))) %>% 
  rename_with(.cols = sort(omniv), ~c(paste0("_omniv",sprintf("%03d", 1:length(omniv))))) %>% 
  rename_with(.cols = 1:1000, ~c(paste0(sprintf("%04d", 1:1000), .)))
rownames(L) = colnames(L)

L = as.matrix(L)

heatweb(L)



# create the 0/1 version of the food web
fw <- L
fw[fw > 0] <- 1

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
      if(assembly:::.components(c(local_producers,local_consumers), fw) == 1 &
         length(assembly:::.find_isolated(c(local_producers,local_consumers), fw)) == 0
      ) break()
    }
    
    
    # subset the meta-fooweb to get the local foodweb
    local_fw = L[c(local_producers, local_consumers),
                 c(local_producers, local_consumers)]
    # add it to the list
    local_fws[[length(local_fws) + 1]] = local_fw
    
    
    
    
    #show_fw(vegan::decostand(local_fw,"pa"), title = "L-matrix model food web")
    
    #show_graph(colnames(local_fw), fw)
  }
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(k, '/', 4000))
  #Sys.sleep(.05)
  if (k == 4000) cat('- Done!')
}



late_succession = vector(mode = "list", length = length(local_fws))
for (m in 1:length(late_succession)) { 
  
  sp_sim <- bride_of_similarity_filtering(colnames(local_fws[[m]]), 
                                          fw, 
                                          t = .1, 
                                          max.iter = 500) %>% 
    sort()
  
  late_succession[[m]] <- L[sp_sim,
                            sp_sim]
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m, '/', length(late_succession)))
  #Sys.sleep(.05)
  if (m == length(local_fws)) cat('- Done!')
  
}



show_fw(vegan::decostand(local_fws[[1000]],"pa"), title = "L-matrix model food web")
show_fw(vegan::decostand(late_succession[[1000]],"pa"), title = "L-matrix model food web")



# saveRDS(late_succession, file="late_succession20220502.RData")
# saveRDS(local_fws, file="local_fws20220502.RData")

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


bride_of_similarity_filtering <- function (sp.names, metaweb, t = 0, 
                                           method = "jaccard", stat = "mean", 
                                           max.iter = 1000) {
  if (t == 0) 
    message("Temperature 't' = 0; this is a purely deterministic filtering")
  isolated <- assembly:::.find_isolated(sp.names, metaweb)
  if (length(isolated) > 0)
    stop("Isolated species detected in input")
  new_sp <- sp.names
  for (i in seq_len(max.iter)) new_sp <- moov(new_sp, metaweb, 
                                              t, method, stat)
  if (length(sp.names) != length(new_sp)) {
    stop("Number of species changed")
  }
  isolated <- assembly:::.find_isolated(new_sp, metaweb)
  print(isolated)
  # if (length(isolated) > 0) 
  #   stop("Isolated species detected in output")
  if (assembly:::.components(new_sp, metaweb) > 1) 
    warning("Isolated component detected in output")
  return(new_sp)
}


moov <- function (sp.names, metaweb, t = 0, method = "jaccard", stat = "mean") 
{
  g <- graph_from_adjacency_matrix(metaweb[sp.names, sp.names])
  consumers <- intersect(sp.names, assembly:::.consumers(metaweb))
  simil <- similarity(g, vids = which(sp.names %in% consumers), 
                      method = method, mode = "all")
  diag(simil) <- NA
  prob_removed <- apply(simil, MARGIN = 2, stat, na.rm = TRUE)
  
  while(TRUE) { 
    remove <- sample(consumers, size = 1, prob = prob_removed)
    repl <- assembly:::.find_replacements(sp.names, remove, metaweb, keep.n.basal = TRUE)
    new.sp <- union(setdiff(sp.names, remove), repl)
    if(assembly:::.components(new.sp, metaweb)==1 &
       length(assembly:::.find_isolated(new.sp, metaweb)) == 0 ) break()
  }
  #  if (length(.find_isolated(new.sp, metaweb) > 0)) 
  #    return(sp.names)
  new.g <- graph_from_adjacency_matrix(metaweb[new.sp, new.sp])
  consumers <- intersect(new.sp, assembly:::.consumers(metaweb))
  new.simil <- similarity(new.g, vids = which(new.sp %in% 
                                                consumers), method = method, mode = "all")
  diag(new.simil) <- NA
  if (stat == "mean") {
    simil <- mean(simil, na.rm = TRUE)
    new.simil <- mean(new.simil, na.rm = TRUE)
  }
  else if (stat == "sum") {
    simil <- sum(simil, na.rm = TRUE)
    new.simil <- sum(new.simil, na.rm = TRUE)
  }
  else {
    stop("'stat' must be one of c('mean', 'sum')")
  }
  if (assembly::metropolis.hastings(simil, new.simil, t = t)) {
    return(new.sp)
  }
  else {
    return(sp.names)
  }
}
