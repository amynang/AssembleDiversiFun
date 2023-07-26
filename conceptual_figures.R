library(tidyverse) 
#library(tidybayes)
library(igraph) # creating networks
library(ggraph) # plotting networks
library(see) # for the half-violin
library(scico) # for the palette
#library(patchwork)
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
# set.seed(321)
# 
# # create dummy data
# plant.niche <- data.frame(Succession = rep(c("Early","Late"),each=1e4),
#                           species = rep(c("sp. 1", "sp. 2"),5e3))
# 
# # draw values from beta to simulate large and small niche overlap
# plant.niche$niche[plant.niche$Succession=="Early" & 
#                     plant.niche$species=="sp. 1"] = rbeta(5e3, 10,10) +.1
# plant.niche$niche[plant.niche$Succession=="Early" & 
#                     plant.niche$species=="sp. 2"] = rbeta(5e3, 10,10) -.1
# plant.niche$niche[plant.niche$Succession=="Late" & 
#                     plant.niche$species=="sp. 1"] = rbeta(5e3, 8,2)
# plant.niche$niche[plant.niche$Succession=="Late" & 
#                     plant.niche$species=="sp. 2"] = rbeta(5e3, 2,8)
# 
# # choose colours
# scico(20, palette = 'bamako')
# scales::show_col(scico(20, palette = 'bamako'))
# 
# pn = ggplot(plant.niche, aes(x = Succession, y=niche, fill = species)) +
#   geom_violinhalf() +
#   geom_hline(yintercept=c(.25, .75), linetype='dashed', color="grey") +
#   theme_modern() +
#   labs(title = ~bold("Niche differentiation")) +
#   theme(  plot.title = element_text(size = 12, hjust = 0.5),
#         axis.title.x = element_text(size = 10),
#         axis.title.y = element_text(size = 10),
#          axis.text.y = element_text(size = 10)) +
#   labs(y = "Niche gradient") +
#   scale_fill_manual(values=c("#36622C","#DAC051"))# + coord_fixed()

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
comp[comp$scenario == "shift" & comp$species=="Sp.1", ]$niche =    rnorm(n, 0,1.5) + 3
comp[comp$scenario == "shift" & comp$species=="Sp.2", ]$niche =    rnorm(n, 0,1.5) - 3

comp$scenario = factor(comp$scenario, levels = rev(c("overlap","restriction","shift")))

pn = ggplot(comp, aes(x=scenario,
                      y=niche,
                      fill=species))+
  geom_violinhalf(position = "identity",
                  trim = F,
                  alpha = .75) +
  #stat_summary(fun = "mean",
  #             geom = "point",
  #             size = 4,
  #             shape = 21,
  #             aes(fill = species),
  #             color = "black") +
  theme_modern() +
  labs(title = ~bold("Plant niche differentiation")) +
  theme(  plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_blank(),#element_text(size = 10),
          axis.line.y = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(5.5, 6*5.5, 5.5, 5.5), "pt")) +
  labs(y = "Niche gradient", 
       x = "") +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051")) +
  theme_tufte()
theme_modern()$plot.margin

set.seed(321)
vec = sort(rexp(1e3,.2), decreasing = T) + rnorm(1e3,0,1)
link.sim = data.frame(Linkage_Similarity = (vec - min(vec)) / (max(vec) - min(vec)),
                      Iterations = 1:1e3)
ls = ggplot(link.sim,aes(Iterations,Linkage_Similarity)) +
  geom_smooth(span = .3, se = FALSE, color="#8F403D")+
  theme_modern() +
  labs(title = ~bold("Animal niche differentiation")) +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
         axis.text.x = element_text(size = 10),
         axis.text.y = element_text(size = 10)) +
  labs(y = "Linkage Similarity")# + coord_fixed(1e3)

set.seed(321)
overlap = heatweb(competition.N(.5,.6, plants = 8)) + coord_fixed()
set.seed(321)
constriction = heatweb(competition.N(.8, 1, plants = 8)) + coord_fixed()
set.seed(321)
expansion = competition.N(.8, 1, plants = 8)
set.seed(321)
diag(expansion) = diag(competition.N(.5,.6, plants = 8))
expansion = heatweb(expansion) + coord_fixed()


# # interaction matrix of the metacommunity
# data <- matrix(sample(0:1, 2500, replace=TRUE, prob=c(0.8,0.2)), nrow=50)
# network <- graph_from_adjacency_matrix(matrix(sample(0:1, 2500, replace=TRUE, 
#                                                      prob=c(0.8,0.2)), nrow=50), 
#                                        mode='undirected', diag=F )
# #the network
# plot(network, 
#      layout=layout.circle,
#      vertex.label=NA,
#      vertex.size=10,
#      vertex.color=sample(scico(50, palette = 'bilbao', 
#                                begin = .2, 
#                                end = .8),
#                          50, replace = F))
# regional.fw <- function(producers = 8, consumers = 30) {
#   
#   # number of all species
#   species = consumers + producers
#   # # matrix of coordinates
#   # lay<-matrix(nrow=species,ncol=2)
#   # lay[,1]<-c(seq(from = 1, to = 20, by = 20/producers),runif(consumers,1,20))
#   # lay[,2]<-c(rep(0,producers),runif(consumers,1,20))
#   
#   # plot
#   p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(sample(0:1, species*species, replace=TRUE, 
#                                                                         prob=c(0.8,0.2)), nrow=species), 
#                                                           mode='undirected', diag=F ),
#                       layout="circle") + 
#     geom_edge_link(colour = "grey") + 
#     geom_node_point(fill = sample(c(scico::scico(producers, palette = 'bamako',
#                                                  end = .8), 
#                                     scico::scico(consumers, palette = 'bilbao', 
#                                                  begin = .2, 
#                                                  end = .8)),
#                                   species, replace = F),
#                     colour = "black",
#                     size = 5,
#                     shape = 21,
#                     stroke = 1) + 
#     theme(legend.position = 'none') +
#     theme_graph(background = "white")
#   return(p)
#   
# }





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

# similar links
module1 <- function(producers = 2, consumers = 2) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(1,2,2,1)
  lay[,2]<-c(1,1,2,2)
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(c(0,0,1,1,
                                                                   0,0,1,1,
                                                                   0,0,0,0,
                                                                   0,0,0,0),
                                                                 ncol = 4),
                                                          mode='undirected', diag=F),
                      layout=lay)+
    geom_edge_link(colour = "black",
                   width = 2) +
    geom_node_point(fill = c("#DAC051","#36622C",
                             "#954B49","#B1A273"),
                    colour = "black",
                    size = 6,
                    shape = 21,
                    stroke = 2) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}

# dissimilar links
module2 <- function(producers = 2, consumers = 2) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(1,2,2,1)
  lay[,2]<-c(1,1,2,2)
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(c(0,0,0,1,
                                                                   0,0,1,0,
                                                                   0,0,0,0,
                                                                   0,0,0,0),
                                                                 ncol = 4),
                                                          mode='undirected', diag=F),
                      layout=lay)+
    geom_edge_link(colour = "black",
                   width = 2) +
    geom_node_point(fill = c("#DAC051","#36622C",
                             "#954B49","#B1A273"),
                    colour = "black",
                    size = 6,
                    shape = 21,
                    stroke = 2) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}


local.fw(250,750)

set.seed(321)
fw0 = plot.regional(100)
fw1 = local.fw( 2,30)
fw2 = local.fw( 2,30)
fw3 = local.fw( 4,30)
fw4 = local.fw( 4,30)
fw5 = local.fw( 8,30)
fw6 = local.fw( 8,30)
fw7 = local.fw(16,30)
fw8 = local.fw(16,30)


early <- plot_grid(fw1,fw3,fw5,fw7, ncol = 1,
                   labels = c("b","","",""),
                   label_size = 15)  + draw_plot_label( "Local communities", size = 12, x = .06)
#late  <- plot_grid(fw2,fw4,fw6,fw8, ncol = 1)  + draw_plot_label("Late succession", size = 12, x = .12)
regional <- plot_grid(fw0, ncol = 1,
                      labels = "a",
                      label_size = 15)     + draw_plot_label("Regional pool", size = 12, x = .33)
#comp <- plot_grid(a.early, a.late, ncol = 2) #+ draw_plot_label("Producer competition", size = 12, x = .16)
mods <- plot_grid(NULL, module1(2,2),module2(2,2), ncol = 3,
                  rel_widths = c(.2,1,1))
sim <- plot_grid(pn, NULL, ls, mods, ncol = 1,
                 labels = c("c","","d",""),
                 label_size = 15,
                 rel_heights = c(2,.5, 1.5,1)) #+ draw_plot_label("Niche differentiation", size = 12, x = .16)
#(concept <- plot_grid(regional, early, sim, late, ncol = 4))
(concept <- plot_grid(regional, early, sim, ncol = 3, rel_widths = c(2,1,1)))

# save the plot
ggsave("Concept2.png",
       bg="white",
       scale = 2,
       width = 173,
       height = 173/1.5,
       units = "mm",
       dpi = 300)


ggsave("Concept.png", concept, bg="white",
       width = 15,
       height = 7.5,
       units = "in",
       dpi = 300)





ggraph(igraph::graph_from_adjacency_matrix(matrix(sample(0:1, 100, replace=TRUE, 
                                                         prob=c(0.8,0.2)), nrow=10), 
                                           mode='undirected', diag=F)) + 
  geom_edge_link(colour = "grey") + 
  geom_node_point(fill = "red", 
                  colour = "black", 
                  size = 10) + 
  theme(legend.position = 'none') +
  theme_graph(background = "white")
