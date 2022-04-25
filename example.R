library(tidyverse)
library(assembly)
library(ATNr)
set.seed(123)

n_species <- 1000
n_basal <- 100
# body mass of species
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
# 2) create the food web
# create the L matrix
L <- create_Lmatrix(masses, 
                    n_basal, 
                    Ropt = 10, #100, #3.98
                    gamma = 2, 
                    th = 0.01)
diag(L) = 0
colnames(L) = c(paste0("broducer",sprintf("%04d", 1:n_basal)),
                paste0("consumer",sprintf("%04d", 1:(n_species-n_basal)))
)

rownames(L) = colnames(L)

# create the 0/1 version of the food web
fw <- L
fw[fw > 0] <- 1


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

show_fw(fw, title = "L-matrix model food web")

producer_diversity = c(2,4,8,16)

local_fws = vector(mode = "list")

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
  local_fw = fw[c(local_producers, local_consumers),
                c(local_producers, local_consumers)]
  
  # add that foodweb to the list
  local_fws[[length(local_fws) + 1]] = local_fw
  
  show_fw(local_fw, title = "L-matrix model food web")
  
  #show_graph(colnames(local_fw), fw)
}

# component doublecheck
gg = vector(mode = "numeric", length=length(local_fws))
for (g in 1:length(local_fws)) {
  gg[g] = assembly:::.components(colnames(local_fws[[g]]), fw)
}
length(which(gg>1))


late_succession = vector(mode = "list", length = length(local_fws))

for (m in 1:length(late_succession)) { 
  
  sp_sim = similarity_filtering(colnames(local_fws[[m]]), 
                                fw, 
                                t = 0, 
                                max.iter = 1000) %>% 
    sort()
  
  late_succession[[m]] = fw[sp_sim,
                            sp_sim]
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m, '/160'))
  #Sys.sleep(.05)
  if (m == length(local_fws)) cat('- Done!')
  
}


# seems to select "clones" of existing species
show_fw(local_fws[[1]], title = "L-matrix model food web")
show_fw(late_succession[[1]], title = "L-matrix model food web")


# average pairwise similarity goes up
similarity.jaccard(graph_from_adjacency_matrix(local_fws[[1]][3:62,3:62])) %>% 
  na_if(., 0) %>% 
  mean(na.rm = T)
similarity.jaccard(graph_from_adjacency_matrix(late_succession[[1]][3:62,3:62])) %>% 
  na_if(., 0) %>% 
  mean(na.rm = T)

# the sum of pairwise similarities goes down
similarity.jaccard(graph_from_adjacency_matrix(local_fws[[1]][3:62,3:62])) %>% 
  sum()
similarity.jaccard(graph_from_adjacency_matrix(late_succession[[1]][3:62,3:62])) %>% 
  sum()
