library(tidyverse)
library(assembly)
library(ATNr)

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
                    Ropt = 100, #3.98
                    gamma = 2, 
                    th = 0.01)

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
consumer_diversity = c(30,40,50,60)
# an empty list 
local_fws = vector(mode = "list")

# k iterations
for (k in 1:10) { 
  for (j in consumer_diversity) {
    
    for (i in producer_diversity) {
      
      while(TRUE){ # while loop skips communities with isolated components
        # sample n=j local consumers from the regional pool
        local_consumers = sample(colnames(fw)[(n_basal+1):n_species], j) %>% 
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
      
      # Do I need this?
      # try(
      #sp_resource <- resource_filtering(local_fw, fw, keep.n.basal = TRUE)
      #   )
      # 
      # local_fw = fw[sp_resource,
      #               sp_resource]
      
      # add that foodweb to the list
      local_fws[[length(local_fws) + 1]] = local_fw
      
      show_fw(local_fw, title = "L-matrix model food web")
      
      #show_graph(colnames(local_fw), fw)
    }
  }
}
gg = vector(mode = "numeric", length=length(local_fws))
for (g in 1:length(local_fws)) {
  gg[g] = assembly:::.components(colnames(local_fws[[g]]), fw)
}
length(which(gg>1))
# loc_fws = local_fws[-which(gg>1)]
# View(local_fws[-which(gg>1)])

late_succession = vector(mode = "list", length = length(local_fws))

for (m in 1:length(late_succession)) { 
  
  sp_sim <- bride_of_similarity_filtering(colnames(local_fws[[m]]), 
                                          fw, 
                                          t = 0, 
                                          max.iter = 100) %>% 
    sort()
  
  late_succession[[m]] <- fw[sp_sim,
                             sp_sim]
  
  
  cat('\014')
  #cat(paste0(round((m/1600)*100), '%'))
  cat(paste0(m, '/160'))
  #Sys.sleep(.05)
  if (m == length(local_fws)) cat('- Done!')
  
}
beepr::beep(9)


show_fw(local_fws[[5]], title = "L-matrix model food web")
show_fw(late_succession[[5]], title = "L-matrix model food web")




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
