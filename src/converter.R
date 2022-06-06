library(igraph)
library(osmar)

setwd("C:/Users/olive/eclipse-workspace/MaxGapMinimization/src")

map = get_osm(complete_file(), source = osmsource_file("map.osm"))

for (i in range(0, map.nodes.size -1)) {
  print(map.nodes[i])
}

plot(map)
graph = as_igraph(map)
