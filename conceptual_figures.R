library(tidyverse)
library(tidybayes)
library(igraph)
library(ggraph)
library(see) # for the half-violin
library(scico) # for the palette
library(patchwork)


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

ggplot(plant.niche, aes(x = Succession, y=niche, fill = species)) +
  geom_violinhalf() +
  geom_hline(yintercept=c(.25, .75), linetype='dashed') +
  theme_modern() +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)) +
  labs(y = "Niche gradient") +
  scale_fill_manual(values=c("#36622C","#DAC051"))



# interaction matrix
data <- matrix(sample(0:1, 2500, replace=TRUE, prob=c(0.8,0.2)), nrow=50)
network <- graph_from_adjacency_matrix(matrix(sample(0:1, 2500, replace=TRUE, 
                                                     prob=c(0.8,0.2)), nrow=50), 
                                       mode='undirected', diag=F )
plot(network, 
     layout=layout.circle,
     vertex.label=NA,
     vertex.size=10,
     vertex.color=sample(scico(50, palette = 'bilbao', 
                               begin = .2, 
                               end = .8),
                         50, replace = F))



create_Lmatrix(c(10^runif(8,-9,-3),
                 10^runif(30,-9,3)), 
               8, Ropt = 3.98, gamma = 2, th = 0.01)

consumers = 30
producers = 8
species = consumers + producers

lay<-matrix(nrow=species,ncol=2)
lay[,1]<-c(seq(from = 1, to = 10, by = 10/producers),runif(consumers,1,10))
lay[,2]<-c(rep(0,producers),runif(consumers,1,10))

plot(graph_from_adjacency_matrix(vegan::decostand(create_Lmatrix(10 ^ c(sort(runif(producers, -9, -3)),
                                                                        sort(runif(consumers, -9, 3))), 
                                                                 producers, 
                                                                 Ropt = 3.98, 
                                                                 gamma = 2, 
                                                                 th = 0.01),
                                                  "pa"), 
                                 mode='undirected', diag=F),
     layout=lay,
     vertex.label=NA,
     vertex.size=10,
     vertex.color=sample(scico(species, palette = 'bilbao', 
                               begin = .2, 
                               end = .8),
                         species, replace = F))


local.fw <- function(producers = 8, consumers = 30) {
  
  # number of all species
  species = consumers + producers
  # matrix of coordinates
  lay<-matrix(nrow=species,ncol=2)
  lay[,1]<-c(seq(from = 1, to = 10, by = 10/producers),runif(consumers,1,10))
  lay[,2]<-c(rep(0,producers),runif(consumers,1,10))
  
  # plot
  p = plot(igraph::graph_from_adjacency_matrix(vegan::decostand(create_Lmatrix(10 ^ c(sort(runif(producers, -9, -3)),
                                                                                      sort(runif(consumers, -9, 3))), 
                                                                               producers, 
                                                                               Ropt = 3.98, 
                                                                               gamma = 2, 
                                                                               th = 0.01),
                                                                "pa"), 
                                               mode='undirected', diag=F),
           layout=lay,
           vertex.label=NA,
           vertex.size=10,
           vertex.color=c(sample(scico::scico(producers, palette = 'bamako'),
                                 producers, replace = F),
                          sample(scico::scico(consumers, palette = 'bilbao', 
                                       begin = .2, 
                                       end = .8),
                                 consumers, replace = F)))
  return(p)
  
}

ggraph(igraph::graph_from_adjacency_matrix(vegan::decostand(create_Lmatrix(10 ^ c(sort(runif(producers, -9, -3)),
                                                                                  sort(runif(consumers, -9, 3))), 
                                                                           producers, 
                                                                           Ropt = 3.98, 
                                                                           gamma = 2, 
                                                                           th = 0.01),
                                                            "pa"), 
                                           mode='undirected', diag=F),
       layout=lay)+ 
  geom_edge_link(colour = "grey") + 
  geom_node_point(colour = c(sample(scico::scico(producers, palette = 'bamako'),
                                    producers, replace = F),
                             sample(scico::scico(consumers, palette = 'bilbao', 
                                                 begin = .2, 
                                                 end = .8),
                                    consumers, replace = F)),
                  size = 10) + 
  theme(legend.position = 'none') +
  theme_graph(background = "white")



fw1 = local.fw(2,30)
fw2 = local.fw(2,30)
fw3 = local.fw(4,30)
fw4 = local.fw(4,30)
fw5 = local.fw(8,30)
fw6 = local.fw(8,30)

par(mfrow=c(3,2))
local.fw(2,30)
local.fw(2,30)
local.fw(4,30)
local.fw(4,30)
local.fw(8,30)
local.fw(8,30)

(fw1 + fw2)/(fw3 + fw4)/(fw5 + fw6)


ggraph(igraph::graph_from_adjacency_matrix(matrix(sample(0:1, 100, replace=TRUE, 
                                                         prob=c(0.8,0.2)), nrow=10), 
                                           mode='undirected', diag=F)) + 
  geom_edge_link(colour = "grey") + 
  geom_node_point(fill = "red", 
                  colour = "black", 
                  size = 10) + 
  theme(legend.position = 'none') +
  theme_graph(background = "white")
