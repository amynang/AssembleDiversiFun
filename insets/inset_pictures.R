library(tidyverse)
#library(tidybayes)
library(igraph) # creating networks
library(ggraph) # plotting networks
library(see) # for the half-violin
library(scico) # for the palette
library(patchwork)


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
                    size = 12,
                    shape = 21,
                    stroke = 2) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}
# module1(2,2)

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
                    size = 12,
                    shape = 21,
                    stroke = 2) +
    theme(legend.position = 'none',
          plot.margin = unit(c(0, 0, -4, 0), "cm")) +
    theme_graph(background = "transparent") + 
    coord_fixed(clip = "off")
  return(p)
  
}
# module2(2,2)

n = 1e4
rlnormtrunc.intuitive = function(n, m, s, p=.9) {
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
comp = data.frame(scenario=rep(c("overlap","constriction","expansion"),
                               each = n*2),
                  species = rep(rep(c("Sp.1","Sp.2"), each = n), 3),
                  niche = NA)
comp[comp$scenario == "overlap" & comp$species == "Sp.1", ]$niche =    rnorm(n, 0,1) + 1
comp[comp$scenario == "overlap" & comp$species == "Sp.2", ]$niche =    rnorm(n, 0,1) - 1
comp[comp$scenario == "constriction" & comp$species=="Sp.1", ]$niche = -rlnormtrunc.intuitive(n, 1.65,1.65) +3
comp[comp$scenario == "constriction" & comp$species=="Sp.2", ]$niche = rlnormtrunc.intuitive(n, 1.65,1.65) -3
comp[comp$scenario == "expansion" & comp$species=="Sp.1", ]$niche =    rnorm(n, 0,1.5) + 3
comp[comp$scenario == "expansion" & comp$species=="Sp.2", ]$niche =    rnorm(n, 0,1.5) - 3

comp$scenario = factor(comp$scenario, levels = rev(c("overlap","constriction","expansion")))

p1 = ggplot(comp[comp$scenario == "overlap",], aes(x=scenario,
                                                   y=niche,
                                                   fill=species))+
  geom_violinhalf(position = position_nudge(x=-.5,y=0),
                  trim = F,
                  alpha = .75) +
  theme_minimal() +
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(-10, 2, -1, 0), "cm")) +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051"))

p2 = ggplot(comp[comp$scenario == "constriction",], aes(x=scenario,
                                                   y=niche,
                                                   fill=species))+
  geom_violinhalf(position = position_nudge(x=-.5,y=0),
                  trim = F,
                  alpha = .75) +
  theme_minimal() +
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(-10, 5, -1, 0), "cm")) +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051"))

p3 = ggplot(comp[comp$scenario == "expansion",], aes(x=scenario,
                                                     y=niche,
                                                     fill=species))+
  geom_violinhalf(position = position_nudge(x=-.5,y=0),
                  trim = F,
                  alpha = .75) +
  theme_minimal() +
  theme(legend.position="none",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(-10, 2, -1, 0), "cm")) +
  coord_flip()+
  scale_fill_manual(values=c("#36622C","#DAC051"))

a1 = module1(2,2)
a2 = module2(2,2)



(p1 + 
    inset_element(a1, 
                    right = 1.8,
                    bottom = 0, 
                    left = 0,
                    top = .6) & 
      plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))

(p1 + 
    inset_element(a2, 
                  right = 1.8,
                  bottom = 0, 
                  left = 0,
                  top = .6) & 
    plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))

(p2 + 
    inset_element(a1, 
                  right = 1.8,
                  bottom = 0, 
                  left = 0,
                  top = .6) & 
    plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))
(p2 + 
    inset_element(a2, 
                  right = 1.8,
                  bottom = 0, 
                  left = 0,
                  top = .6) & 
    plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))

(p3 + 
    inset_element(a1, 
                  right = 1.8,
                  bottom = 0, 
                  left = 0,
                  top = .6) & 
    plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))
(p3 + 
    inset_element(a2, 
                  right = 1.8,
                  bottom = 0, 
                  left = 0,
                  top = .6) & 
    plot_layout(widths = c(10,1)) &
    theme(plot.background = element_rect(fill='transparent', 
                                         color=NA)))

# save each plot
ggsave("ll1.png",
       scale = 3,
       width = 80,
       height = 40,
       units = "mm",
       dpi = 300)
