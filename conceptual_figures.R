library(tidyverse)
library(tidybayes)
library(igraph)
library(ggraph)
library(see) # for the half-violin
library(scico) # for the palette
library(patchwork)
library(cowplot)


set.seed(321)

# create dummy data
plant.niche <- data.frame(Succession = rep(c("Early","Late"),each=1e4),
                          species = rep(c("sp. 1", "sp. 2"),5e3))

# draw values from beta to simulate large and small niche overlap
plant.niche$niche[plant.niche$Succession=="Early" & 
                    plant.niche$species=="sp. 1"] = rbeta(5e3, 10,10) +.1
plant.niche$niche[plant.niche$Succession=="Early" & 
                    plant.niche$species=="sp. 2"] = rbeta(5e3, 10,10) -.1
plant.niche$niche[plant.niche$Succession=="Late" & 
                    plant.niche$species=="sp. 1"] = rbeta(5e3, 8,2)
plant.niche$niche[plant.niche$Succession=="Late" & 
                    plant.niche$species=="sp. 2"] = rbeta(5e3, 2,8)

# choose colours
scico(20, palette = 'bamako')
scales::show_col(scico(20, palette = 'bamako'))

pn = ggplot(plant.niche, aes(x = Succession, y=niche, fill = species)) +
  geom_violinhalf() +
  geom_hline(yintercept=c(.25, .75), linetype='dashed', color="grey") +
  theme_modern() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  labs(y = "Niche gradient") +
  scale_fill_manual(values=c("#36622C","#DAC051"))



set.seed(321)
vec = sort(rexp(1e3,.2), decreasing = T) + rnorm(1e3,0,1)
link.sim = data.frame(Linkage_Similarity = (vec - min(vec)) / (max(vec) - min(vec)),
                      Iterations = 1:1e3)
ls = ggplot(link.sim,aes(Iterations,Linkage_Similarity)) +
  geom_smooth(span = .3, se = FALSE, color="#8F403D")+
  theme_modern() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  labs(y = "Linkage Similarity")


# interaction matrix of the metacommunity
data <- matrix(sample(0:1, 2500, replace=TRUE, prob=c(0.8,0.2)), nrow=50)
network <- graph_from_adjacency_matrix(matrix(sample(0:1, 2500, replace=TRUE, 
                                                     prob=c(0.8,0.2)), nrow=50), 
                                       mode='undirected', diag=F )
#the network
plot(network, 
     layout=layout.circle,
     vertex.label=NA,
     vertex.size=10,
     vertex.color=sample(scico(50, palette = 'bilbao', 
                               begin = .2, 
                               end = .8),
                         50, replace = F))
regional.fw <- function(producers = 8, consumers = 30) {
  
  # number of all species
  species = consumers + producers
  # # matrix of coordinates
  # lay<-matrix(nrow=species,ncol=2)
  # lay[,1]<-c(seq(from = 1, to = 20, by = 20/producers),runif(consumers,1,20))
  # lay[,2]<-c(rep(0,producers),runif(consumers,1,20))
  
  # plot
  p <- ggraph::ggraph(igraph::graph_from_adjacency_matrix(matrix(sample(0:1, species*species, replace=TRUE, 
                                                                        prob=c(0.8,0.2)), nrow=species), 
                                                          mode='undirected', diag=F ),
                      layout="circle") + 
    geom_edge_link(colour = "grey") + 
    geom_node_point(fill = sample(c(scico::scico(producers, palette = 'bamako',
                                                 end = .8), 
                                    scico::scico(consumers, palette = 'bilbao', 
                                                 begin = .2, 
                                                 end = .8)),
                                  species, replace = F),
                    colour = "black",
                    size = 5,
                    shape = 21,
                    stroke = 1) + 
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
    theme_graph(background = "white")
  return(p)
  
}





local.fw(250,750)

set.seed(321)
fw0 = regional.fw(13,37)
fw1 = local.fw( 2,30)
fw2 = local.fw( 2,30)
fw3 = local.fw( 4,30)
fw4 = local.fw( 4,30)
fw5 = local.fw( 8,30)
fw6 = local.fw( 8,30)
fw7 = local.fw(16,30)
fw8 = local.fw(16,30)

(plot_spacer()/fw0) + ((fw1 + fw2)/(fw3 + fw4)/(fw5 + fw6)/(fw7 + fw8))



local <- plot_grid(fw1,fw2,fw3,fw4,fw5,fw6,fw7,fw8, ncol = 2)
regional <- plot_grid(fw0,NULL, ncol = 1)
sim <- plot_grid(pn,ls,NULL,ncol = 1)
plot_grid(regional, local, sim, ncol = 3)



+ plot_annotation(tag_levels = 'a')


ggraph(igraph::graph_from_adjacency_matrix(matrix(sample(0:1, 100, replace=TRUE, 
                                                         prob=c(0.8,0.2)), nrow=10), 
                                           mode='undirected', diag=F)) + 
  geom_edge_link(colour = "grey") + 
  geom_node_point(fill = "red", 
                  colour = "black", 
                  size = 10) + 
  theme(legend.position = 'none') +
  theme_graph(background = "white")
