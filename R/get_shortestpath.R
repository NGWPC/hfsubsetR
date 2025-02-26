# Caleb Novinger 
# Shortest Path Subset - HFSubsetR
# Install and load required packages
install.packages(c("sf", "sfnetworks", "tidygraph", "dplyr","remotes"))
remotes::install_github("owp-spatial/hfsubsetR")
remotes::install_github("NOAA-OWP/hydrofabric")
library(sf)
library(tidyverse)
library(igraph)
library(dplyr)
library(tidygraph)
# Define the function
get_shortest_path <- function(gpkg_path, start_hf_id, end_hf_id, output_gpkg_path) {
  # Read the layers from the GeoPackage
  flowpaths <- st_read(gpkg_path, layer = "flowpaths", quiet = TRUE)
  divides <- st_read(gpkg_path, layer = "divides", quiet = TRUE)
  nexus <- st_read(gpkg_path, layer = "nexus", quiet = TRUE)
  network <- st_read(gpkg_path, layer = "network", quiet = TRUE)
  # Ensure hf_id, id, and toid are character strings in the network layer
  network <- network %>%
    mutate(
      hf_id = as.character(hf_id), 
      id = as.character(id), 
      toid = as.character(toid)
    )
  # Find corresponding "id" values for given hf_id inputs
  start_node <- network %>% filter(hf_id == start_hf_id) %>% pull(id)
  end_node <- network %>% filter(hf_id == end_hf_id) %>% pull(id)
  
  if (length(start_node) == 0) stop(paste("Start hf_id", start_hf_id, "not found in the network!"))
  if (length(end_node) == 0) stop(paste("End hf_id", end_hf_id, "not found in the network!"))
  # Create a directed graph
  graph <- as_tbl_graph(network, directed = TRUE) %>%
    activate("edges") %>%
    mutate(weight = ifelse(is.na(lengthkm), 1, lengthkm))  # Assign weight based on lengthkm
  # Convert graph to igraph
  igraph_graph <- as.igraph(graph)
  # Check if start and end nodes exist in the graph
  all_nodes <- V(igraph_graph)$name  
  if (!(start_node %in% all_nodes)) stop(paste("Start node", start_node, "not found in the graph!"))
  if (!(end_node %in% all_nodes)) stop(paste("End node", end_node, "not found in the graph!"))
  # Find shortest path using Dijkstra’s algorithm
  shortest_path <- igraph::shortest_paths(
    igraph_graph,
    from = which(V(igraph_graph)$name == start_node),
    to = which(V(igraph_graph)$name == end_node),
    weights = E(igraph_graph)$weight,
    mode = "out"
  )
  # Extract path nodes
  path_nodes <- V(igraph_graph)[shortest_path$vpath[[1]]]$name
  print(path_nodes)
  
  # Ensure IDs are character for proper filtering
  flowpaths <- flowpaths %>% mutate(id = as.character(id))
  divides <- divides %>% mutate(id = as.character(id))
  nexus <- nexus %>% mutate(id = as.character(id))
  network <- network %>% mutate(id = as.character(id))
  # Filter flowpaths to match the shortest path
  shortest_flowpaths <- flowpaths %>% filter(id %in% path_nodes)
  # Identify divides, nexus, and network features linked to the shortest path flowpaths
  shortest_divides <- divides %>% filter(id %in% shortest_flowpaths$id)
  shortest_nexus <- nexus %>% filter(id %in% c(shortest_flowpaths$id, shortest_flowpaths$toid)) 
  shortest_network <- network %>% filter(id %in% shortest_flowpaths$id | toid %in% shortest_flowpaths$id)
  # Write only the subsetted layers to the output GeoPackage
  st_write(shortest_flowpaths, output_gpkg_path, layer = "shortest_flowpaths", delete_layer = TRUE)
  st_write(shortest_divides, output_gpkg_path, layer = "shortest_divides", delete_layer = FALSE)
  st_write(shortest_nexus, output_gpkg_path, layer = "shortest_nexus", delete_layer = FALSE)
  st_write(shortest_network, output_gpkg_path, layer = "shortest_network", delete_layer = FALSE)
  return(paste("GeoPackage saved to:", output_gpkg_path))
}

get_shortest_path("C:/Users/cnovi/Downloads/conus_nextgen.gpkg", 
                  "22162876", "22144632", "C:/Users/cnovi/Downloads/conus_nextgenOUTPUT.gpkg")