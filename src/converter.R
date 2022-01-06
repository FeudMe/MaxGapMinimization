library(igraph)
library(osmar)

setwd("C:/Users/olive/eclipse-workspace/MaxGapMinimization/src")

file = get_osm(complete_file(), source = osmsource_file("map.osm"))

graph = as.undirected(as_igraph(file))
plot(graph, vertex.label = NA, vertex.size = 0, vertex.shape = 'none', edge.color = "black", edge.)