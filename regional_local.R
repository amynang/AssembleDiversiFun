library(tidyverse)
library(assembly)
library(ATNr)
set.seed(123)

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


# specify number of plants and animals
dim1 = 250
dim2 = 750

# quick and dirty way for nested plant-herbivore interactions
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

# shuffle herbivores so that later on this pattern will be unrelated to bodymasses
rand <- sample(ncol(drum))
#rand
drum = drum[,rand]

# create the interaction matrix
L = matrix(NA,dim1+dim2,dim1+dim2)
# plants don't eat anything
L[,1:dim1] = 0
# put in the herbivory
L[1:dim1,(dim1+1):(dim1+dim2)] = drum#[dim(drum)[1]:1, ]

n_species <- dim1+dim2
n_basal <- dim1
#n_nut <- 16
# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
# 2) create the food web
# create the L matrix
K <- create_Lmatrix(masses, 
                    n_basal, 
                    Ropt = 3.98, #100, #3.98
                    gamma = 2, 
                    th = 0.01)
# predatory links are allometric
L[(dim1+1):(dim1+dim2), ] = K[(dim1+1):(dim1+dim2), ]

# remove cannibalism
diag(L) = 0

# select random subest of consumers to be herbivores
herbiv = sample((dim1+1):(dim1+dim2), 250)
# convert them to veganism
L[(dim1+1):(dim1+dim2), herbiv] = 0

# select a different random subset of consumers to be predators
predat = sample(setdiff((dim1+1):(dim1+dim2), herbiv), 250)
# make them only eat meat
L[1:dim1, predat] = 0

# the remaining 250 animals will be omnivorous
omniv = setdiff(1:1000,c(1:dim1,herbiv,predat))

show_fw(vegan::decostand(L,"pa"))

# species names are 0000group000 where 0000 is the index of the species (relates to bodymass) group is
# plants,herbivores,omnivores,predators and 000 is the within group index
L = as.data.frame(L) %>%  rename_with(.cols = 1:dim1, ~c(paste0("_plant",sprintf("%03d", 1:250)))) %>% 
                          rename_with(.cols = sort(herbiv), ~c(paste0("_herbiv",sprintf("%03d", 1:250)))) %>% 
                          rename_with(.cols = sort(predat), ~c(paste0("_predat",sprintf("%03d", 1:250)))) %>% 
                          rename_with(.cols = sort(omniv), ~c(paste0("_omniv",sprintf("%03d", 1:250)))) %>% 
                          rename_with(.cols = 1:1000, ~c(paste0(sprintf("%04d", 1:1000), .)))
rownames(L) = colnames(L)

L = as.matrix(L)

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
  
  for (i in producer_diversity*3) {
    
    while(TRUE){ # while loop skips communities with isolated components
      # sample n=j local consumers from the regional pool
      local_consumers = sample(colnames(fw)[(n_basal+1):n_species], 100) %>% 
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


late_succession = vector(mode = "list", length = length(local_fws))

for (m in 1:length(late_succession)) { 
  counter = 0
  while(TRUE) { 
    
    # this will give 0x0 matrices in case of error (isolated species or components)
    sp_sim <- tryCatch(similarity_filtering(colnames(local_fws[[m]]), 
                                            fw, 
                                            t = .01, 
                                            max.iter = 500) %>% 
                         sort(),
                       error = function(e) 0)
    
    late_succession[[m]] <- L[sp_sim,
                              sp_sim]
    counter <- counter+1
    if(dim(late_succession[[m]])[1]!=0 | counter==5
       ) break()
  }
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m, '/', length(late_succession)))
  #Sys.sleep(.05)
  if (m == length(local_fws)) cat('- Done!')
  
}
dims = vector(mode = "numeric", length(late_succession))
for (i in 1:length(late_succession)) {
  dims[i] = dim(late_succession[[i]])[1]
}
table(dims) 


show_fw(vegan::decostand(local_fws[[15]],"pa"), title = "L-matrix model food web")
show_fw(vegan::decostand(late_succession[[15]],"pa"), title = "L-matrix model food web")
