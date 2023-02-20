library(tidyverse) 
#library(tidybayes)
library(igraph) # creating networks
library(ggraph) # plotting networks
library(see) # for the half-violin
library(scico) # for the palette
library(grid)
library(cowplot)
library(ATNr)

source("functions.R")

rlnormtrunc.intuitive = function(n, m, s, p=.95) {
  trnc <- EnvStats::rlnormTrunc(n, 
                                meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                sdlog = sqrt(log(1 + (s^2 / m^2))), 
                                min = qlnorm((1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))), 
                                max = qlnorm(1-(1-p)/2, 
                                             meanlog = log(m^2 / sqrt(s^2 + m^2)), 
                                             sdlog = sqrt(log(1 + (s^2 / m^2)))))
  return(trnc)
}
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
    scale_fill_viridis_c(option = "inferno", direction = -1, begin = .2)
  return(heat)
}

# create example regional matrix
create.regional <- function (species = 100) { 
  set.seed(321)
  
  # specify number of plants and plant consumers (herbivores & omnivores)
  dim1 = species/4
  dim2 = species/2
  
  # quick & dirty way for nested plant-consumer interactions
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
  L = matrix(0,species,species)
  # plants don't eat anything
  L[,1:dim1] = 0
  
  # randomly select 2/3 of animals to be plant consumers (incl. omnivores)
  plant.cons = sample((dim1+1):species, dim2) %>% sort()
  
  # put in the herbivory
  L[1:dim1,plant.cons] = drum#[dim(drum)[1]:1, ]
  
  n_species <- species
  n_basal <- dim1
  
  # body mass of species
  # masses <- 10 ^ c(sort(runif(n_basal, 0, 6)),
  #                  sort(runif(n_species - n_basal, 2, 12)))
  
  # # from 1ngram to 1kgram
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
    L[(dim1+1):n_species, i] = L[(dim1+1):n_species, i]*(1-rbernoulli(1, sum(vegan::decostand(L[1:dim1,i],"pa"))/(dim2/2)))
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
  colSums(vegan::decostand(L[,(dim1+1):species],"pa"))
  min(colSums(vegan::decostand(L[,(dim1+1):species],"pa")))
  
  # species names are 0000group000 where 0000 is the index of the species (relates to bodymass) group is
  # plants,herbivores,omnivores,predators and 000 is the within group index
  L = as.data.frame(L) %>%  rename_with(.cols = 1:dim1, ~c(paste0("_plant",sprintf("%03d", 1:dim1)))) %>% 
    rename_with(.cols = sort(herbiv), ~c(paste0("_herbiv",sprintf("%03d", 1:length(herbiv))))) %>% 
    rename_with(.cols = sort(predat), ~c(paste0("_predat",sprintf("%03d", 1:length(predat))))) %>% 
    rename_with(.cols = sort(omniv), ~c(paste0("_omniv",sprintf("%03d", 1:length(omniv))))) %>% 
    rename_with(.cols = 1:species, ~c(paste0(sprintf("%04d", 1:species), .)))
  rownames(L) = colnames(L)
  
  L = as.matrix(L)
  
  return(L)
}


plot.regional <- function(species = 100) { 
  # create example of metacommunity (1/10 the size)
  L = create.regional(species)
  # turn into a graph
  g = graph_from_adjacency_matrix(vegan::decostand(L,"pa"), 
                                  mode='directed', diag=F )
  # turn to edgelist, name columns "from", "to", add type of consumption
  edge = as.data.frame(as_edgelist(g))
  colnames(edge) = c("from","to")
  edge$type = as.factor(ifelse(grepl("pla",edge$from)==TRUE,"herbivory","predation"))
  
  # turn to graph again, now it has "type" as attribute
  g=graph_from_data_frame(edge, directed = TRUE)
  # create layout
  layout <- create_layout(g, layout = 'linear', circular = TRUE)
  
  # generate plot
  p = ggraph::ggraph(layout) + 
    # edges are coloured based on "herbivory"/"predation"
    geom_edge_arc(aes(colour = type, alpha=.5), show.legend = F) +
    # nodes get green hues if plants, sepia for animals
    geom_node_point(fill = c(sample(scico::scico(species/4, palette = 'bamako',
                                                 end = .8),
                                    species/4, replace = F),
                             sample(scico::scico(3*species/4, palette = 'bilbao',
                                                 begin = .2,
                                                 end = .8),
                                    3*species/4, replace = F)),
                    colour = "black",
                    size = 4,
                    shape = 21,
                    stroke = .5) +
    scale_edge_colour_manual(values=c("#7F7659","#B99B76")) +
    coord_fixed(clip = "off") + 
    theme(legend.position = 'none') +
    theme_graph(background = "white")
  
  return(p)
}


local.fw <- function(producers = 8, consumers = 30) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(seq(from = 1, to = 20, by = 20/producers),runif(consumers,1,20))
  lay[,2]<-c(rep(0,producers),runif(consumers,1,20))
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(vegan::decostand(create_Lmatrix(10 ^ c(sort(runif(producers, -9, -3)),
                                                                                                 sort(runif(consumers, -9, 3))), 
                                                                                          producers, 
                                                                                          Ropt = 3.98, 
                                                                                          gamma = 2, 
                                                                                          th = 0.01),
                                                                           "pa"), 
                                                          mode='undirected', diag=F),
                      layout=lay)+ 
    geom_edge_link(colour = "grey") + 
    geom_node_point(fill = c(sample(scico::scico(producers, palette = 'bamako',
                                                 end = .8),
                                    producers, replace = F),
                             sample(scico::scico(consumers, palette = 'bilbao', 
                                                 begin = .2, 
                                                 end = .8),
                                    consumers, replace = F)),
                    colour = "black",
                    size = 3,
                    shape = 21,
                    stroke = .5) + 
    theme(legend.position = 'none') +
    theme_graph(background = "white") + 
    coord_fixed(clip = "off")
  return(p)
  
}

set.seed(321)
n = 1e4

comp = data.frame(scenario=rep(c("overlap","restriction","shift"),
                               each = n*2),
                  species = rep(rep(c("Sp.1","Sp.2"), each = n), 3),
                  niche = NA)
comp[comp$scenario == "overlap" & comp$species == "Sp.1", ]$niche =    rnorm(n, 0,1) + 1
comp[comp$scenario == "overlap" & comp$species == "Sp.2", ]$niche =    rnorm(n, 0,1) - 1
comp[comp$scenario == "restriction" & comp$species=="Sp.1", ]$niche = -rlnormtrunc.intuitive(n, 1.65,1.65) +3
comp[comp$scenario == "restriction" & comp$species=="Sp.2", ]$niche = rlnormtrunc.intuitive(n, 1.65,1.65) -3
comp[comp$scenario == "shift" & comp$species=="Sp.1", ]$niche =    rnorm(n, 0,1) + 2
comp[comp$scenario == "shift" & comp$species=="Sp.2", ]$niche =    rnorm(n, 0,1) - 2

comp$scenario = factor(comp$scenario, levels = rev(c("overlap","restriction","shift")))

comp$scenario = factor(comp$scenario, levels = rev(c("shift","overlap","restriction")))


p1 = ggplot(comp %>% mutate(niche = case_when(scenario == "restriction" ~ as.numeric(NA),
                                              scenario == "shift" ~ as.numeric(NA),
                                              TRUE ~ niche),
                            scenario = fct_relevel(scenario, c("restriction","overlap","shift")))  
            
            , aes(x=scenario,
                  y=niche,
                  fill=species))+
  geom_violinhalf(position = "identity",
                  trim = F,
                  size = 1,
                  alpha = .75) +
  #stat_summary(fun = "mean",
  #             geom = "point",
  #             size = 4,
  #             shape = 21,
  #             aes(fill = species),
  #             color = "black") +
  theme_modern() +
  #labs(title = ~bold("Plant niche differentiation")) +
  theme(  plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(),#element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_blank(),#element_text(size = 10),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(5.5, 6*5.5, 5.5, 5.5), "pt")) +
  labs(y = "Niche gradient", 
       x = "") +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051"))



p2 = ggplot(comp, aes(x=scenario,
                      y=niche,
                      fill=species))+
  geom_violinhalf(position = "identity",
                  trim = F,
                  size = 1,
                  alpha = .75) +
  #stat_summary(fun = "mean",
  #             geom = "point",
  #             size = 4,
  #             shape = 21,
  #             aes(fill = species),
  #             color = "black") +
  theme_modern() +
  #labs(title = ~bold("Plant niche differentiation")) +
  theme(  plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(),#element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_blank(),#element_text(size = 10),
          axis.text.x = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(5.5, 6*5.5, 5.5, 5.5), "pt")) +
  labs(y = "Niche gradient", 
       x = "") +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051"))


l = plot_grid(p2, p2, ncol = 2,
              rel_widths = c(1.5,1.5),
              rel_heights = c(1,2)) +
  annotate("text", x=.75, y=.78, 
           label = "spread out", 
           color = "black", 
           angle = 0, 
           #hjust = 1.5, 
           size = 10, 
           fontface = "plain") +
  annotate("text", x=.75, y=.48, 
           label = "overlapping", 
           color = "black", 
           angle = 0, 
           #hjust = 1.5, 
           size = 10, 
           fontface = "plain") +
  annotate("text", x=.75, y=.178, 
           label = "clumped", 
           color = "black", 
           angle = 0, 
           #hjust = 1.5, 
           size = 10, 
           fontface = "plain") +
  annotate(geom = "curve",
           size = 1,
           x = .45, xend = .525,
           y = .7, yend = .8,
           curvature = -.1,
           arrow = arrow()) +
  annotate(geom = "curve",
           size = 1,
           x = .45, xend = .525,
           y = .5, yend = .5,
           curvature = -.01,
           arrow = arrow()) +
  annotate(geom = "curve",
           size = 1,
           x = .45, xend = .525,
           y = .3, yend = .2,
           curvature = .1,
           arrow = arrow()) +
  annotate(geom = "curve",
           size = 2,
           x = .1, xend = .9,
           y = .1, yend = .1,
           curvature = 0,
           arrow = arrow()) +
  annotate("text", x=.48,y=.79, 
           label = "adaptation", 
           color = "darkgrey", 
           angle = 45, 
           #hjust = 1.5, 
           size = 8, 
           fontface = "bold") +
  annotate("text", x=.48,y=.55, 
           label = "null", 
           color = "darkgrey", 
           angle = 0, 
           #hjust = 1.5, 
           size = 8, 
           fontface = "bold") +
  annotate("text", x=.48,y=.35, 
           label = "concentration", 
           color = "darkgrey", 
           angle = 315, 
           #hjust = 1.5, 
           size = 8, 
           fontface = "bold") +
  annotate("text", x=.85,y=.05, 
           label = "time", 
           color = "black", 
           angle = 0, 
           #hjust = 1.5, 
           size = 10, 
           fontface = "plain")


module1 <- function(producers = 2, consumers = 4) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(1,2,.5,1.5,2.5,1.25)
  lay[,2]<-c(1,1,2,2,2,3)
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(c(0,0,1,1,1,0,
                                                                   0,0,0,1,1,0,
                                                                   0,0,0,0,0,1,
                                                                   0,0,0,0,0,1,
                                                                   0,0,0,0,0,1,
                                                                   0,0,0,0,0,0),
                                                                 ncol = 6),
                                                          mode='undirected', diag=F)
                      ,layout=lay
  )+
    geom_edge_link(colour = "black",
                   width = 2) +
    geom_node_point(fill = c("#b3a01f","#36622C",
                                      "#a88464","#99514e",
                                      "#99514e","#6d2120"),
                                      colour = "black",
                    size = 8,
                    shape = 21,
                    stroke = 1) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}

module1(2,4)


module2 <- function(producers = 2, consumers = 4) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(1,2,.5,1.5,1.25)
  lay[,2]<-c(1,1,2,2,3)
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(c(0,0,1,1,0,
                                                                   0,0,0,1,0,
                                                                   0,0,0,0,1,
                                                                   0,0,0,0,1,
                                                                   0,0,0,0,0),
                                                                 ncol = 5),
                                                          mode='undirected', diag=F)
                      ,layout=lay
  )+
    geom_edge_link(colour = "black",
                   width = 2) +
    geom_node_point(fill = c("#b3a01f","#36622C",
                                      "#a88464","#99514e",
                                      "#6d2120"),
                                      colour = "black",
                    size = 8,
                    shape = 21,
                    stroke = 1) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}
module2(2,3)


module3 <- function(producers = 2, consumers = 4) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(1,2,.5,1.5,1.25,2)
  lay[,2]<-c(1,1,2,2,3,3.5)
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(c(0,0,1,1,0,0,
                                                                   0,0,0,1,0,0,
                                                                   0,0,0,0,1,0,
                                                                   0,0,0,0,1,0,
                                                                   0,0,0,0,0,1,
                                                                   0,0,0,0,0,0),
                                                                 ncol = 6),
                                                          mode='undirected', diag=F)
                      ,layout=lay
  )+
    geom_edge_link(colour = "black",
                   width = 2) +
    geom_node_point(fill = c("#b3a01f","#36622C",
                                      "#a88464","#99514e",
                                      "#6d2120","#bab18e"),
                                      colour = "black",
                    size = 8,
                    shape = 21,
                    stroke = 1) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}
module3(2,4)


scico(20, palette = 'bilbao')
scales::show_col(scico(20, palette = 'bilbao'))


reg = plot.regional(100)

void = ggplot() + theme_void()
(r = ggdraw(void + 
              draw_plot(reg, .25, .575, .5, .5) + 
              draw_plot(module1(2,4),   0, .25, .3, .3) +
              draw_plot(module2(2,3), .30, .35, .3, .3) +
              draw_plot(module3(2,4), .60, .45, .3, .3) +
              draw_plot(module1(2,4), .65, .10, .3, .3)) +
    annotate(geom = "curve", # random assembly
             size = 1,
             x = .35, xend = .175,
             y = .75, yend = .55,
             curvature = .4,
             arrow = arrow(type = "closed")) +
    annotate(geom = "curve", # no change
             size = 1,
             x = .3, xend = .65,
             y = .375, yend = .30,
             curvature = 0,
             arrow = arrow()) +
    annotate(geom = "curve", # exclusion
             size = 1,
             x = .3, xend = .375,
             y = .425, yend = .45,
             curvature = 0,
             arrow = arrow()) +
    annotate(geom = "curve", # the one without a name
             size = 1,
             x = .52, xend = .65,
             y = .5, yend = .55,
             curvature = 0,
             arrow = arrow()) +
    annotate(geom = "curve", # colonisation
             size = 1,
             x = .64, xend = .77,
             y = .77, yend = .68,
             curvature = -.2,
             arrow = arrow(type = "closed"))+
    annotate(geom = "curve", # time
             size = 2,
             x = .1, xend = .9,
             y = .1, yend = .1,
             curvature = 0,
             arrow = arrow()) +
    annotate("text", x=.2,y=.75, 
             label = "random \n assembly", 
             color = "darkgrey", 
             angle = 42, 
             #hjust = 1.5, 
             size = 8, 
             fontface = "bold") +
    annotate("text", x=.45,y=.325, 
             label = "no change", 
             color = "darkgrey", 
             angle = 350, 
             #hjust = 1.5, 
             size = 8, 
             fontface = "bold") +
    annotate("text", x=.33,y=.46, 
             label = "exclusion", 
             color = "darkgrey", 
             angle = 16, 
             #hjust = 1.5, 
             size = 8, 
             fontface = "bold") +
    annotate("text", x=.725,y=.775, 
             label = "colonisation", 
             color = "darkgrey", 
             angle = 330, 
             #hjust = 1.5, 
             size = 8, 
             fontface = "bold") +
    annotate("text", x=.75, y=.478, 
             label = "high complementarity", 
             color = "black", 
             angle = 0, 
             #hjust = 1.5, 
             size = 10, 
             fontface = "plain") +
    annotate("text", x=.75, y=.15, 
             label = "low complementarity", 
             color = "black", 
             angle = 0, 
             #hjust = 1.5, 
             size = 10, 
             fontface = "plain") +
    annotate(geom = "point", 
             x=.2478,y=.4109, 
             shape = 1,
             stroke = 1,
             color = "grey", 
             size = 8) +
    annotate("text", x=.85,y=.05, 
             label = "time", 
             color = "black", 
             angle = 0, 
             #hjust = 1.5, 
             size = 10, 
             fontface = "plain"))

lr = plot_grid(l,r,ncol=2,
               labels = c("a","b"),
               label_size = 30)

set.seed(321)
b = plot_grid(nrow = 1,
              local.fw(2,30),
              local.fw(4,30),
              local.fw(6,30),
              local.fw(8,30),
              local.fw(10,30),
              local.fw(12,30),
              local.fw(14,30),
              local.fw(16,30))

plot_grid(lr,b,NULL,
          ncol = 1,
          labels = c("","c",""),
          label_size = 30,
          rel_heights = c(1,.25,.075)) +
  draw_grob(rectGrob(x = .12, y = .426,
                     width = .11,
                     height = .14,
                     gp = gpar(col = "white",
                               fill = alpha("white", 1),
                               alpha = 1,
                               lwd = 3))) +
  draw_grob(rectGrob(x = .125, y = .9,
                     width = .2,
                     height = .14,
                     gp = gpar(col = "white",
                               fill = alpha("white", 1),
                               alpha = 1,
                               lwd = 3))) +
  annotate(geom = "curve", # time
           size = 2,
           x = .05, xend = .95,
           y = .07, yend = .07,
           curvature = 0,
           arrow = arrow()) +
  annotate("text", x=.9,y=.04, 
           label = "plant richness", 
           color = "black", 
           angle = 0, 
           #hjust = 1.5, 
           size = 10, 
           fontface = "plain")

ggsave("Figure1.png",
       bg="white",
       scale = 4,
       width = 173,
       height = 173/2,
       units = "mm",
       dpi = 300)

