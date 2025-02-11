# Caleb Novinger 
# Shortest Path Subset - HFSubsetR
# Create a new R script in the R/ directory (e.g., R/get_shortest_path_subset.R)

install.packages("remotes")
remotes::install_github("owp-spatial/hfsubsetR") # Load necessary libraries
install.packages("dplyr")
install.packages("igraph")
library(sf)  # For spatial data handling
library(dplyr)  # For data manipulation
library(igraph)
library(hfsubsetR)
# Function to calculate the shortest path subset
get_shortest_path_subset <- function(gpkg, layer, start_id, end_id, outfile) {
  
  # Read in the hydrofabric data from the specified layer
  hydrofabric <- sf::st_read(gpkg, layer = network)
  
  # Retrieve the start and end features based on identifiers (as integers)
  start_point <- hydrofabric %>% filter(id == as.integer(start_id))
  end_point <- hydrofabric %>% filter(id == as.integer(end_id))
  
  # Construct a spatial network
  network <- sf::st_as_sf(hydrofabric)
  
  # Use a graph library to calculate the shortest path
  g <- igraph::graph_from_data_frame(as.data.frame(network), directed = FALSE)
  
  # Calculate the shortest path
  path <- igraph::shortest_paths(g, from = as.integer(start_id), to = as.integer(end_id))
  
  # Convert the path to a subset of the hydrofabric
  shortest_path_subset <- hydrofabric[unlist(path$vpath), ]
  
  # Save the resulting subset to a new geopackage
  sf::st_write(shortest_path_subset, outfile)
  
  return(shortest_path_subset)
}

# Example usage:
# get_shortest_path_subset(gpkg = '~/hydrofabric/v2.2/conus_nextgen.gpkg', layer = 'network', start_id = 22162876, end_id = 22144632, outfile = './shortest_path_subset.gpkg')