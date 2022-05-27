
# to plot the trophic matrix
heatweb <- function(mat) {
  heat <- mat %>% #na_if(., 0) %>% 
    as.data.frame() %>%
    rownames_to_column("id") %>%
    pivot_longer(-c(id), names_to = "species", values_to = "strength") %>%
    mutate(species= fct_relevel(species,colnames(mat))) %>%
    ggplot(aes(x=species, y=ordered(id, levels = rev(unique(id))), fill=strength)) + 
    geom_raster() +
    theme_void() +
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

# to create communities of limited linkage similarity
# edited from https://github.com/emilio-berti/assembly/blob/master/R/limiting_similarity_filtering.R
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

# internal function of above, instead of throwing error uses a while loop to 
# ensure no isolated species / components
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


# to generate the plant competition matrix
competition <- function(lower=.8, upper=1, plants) { 
  # create competition matrix
  alpha = matrix(NA, plants, plants)
  # draw diagonals from (lower, upper)
  diag(alpha) = runif(plants, lower,upper)
  # draw off diagonals so that columns sum to 1
  for (i in 1:plants) {
    # replace NAs in each column with values that sum to the complement of that diagonal
    alpha[which(is.na(alpha[,i])),i] = (1-diag(alpha)[i])*brms::rdirichlet(1,rep(2,(plants-1)))
  }
  return(alpha)
}