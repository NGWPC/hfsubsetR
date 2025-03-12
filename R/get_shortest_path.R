#' @title Find the shortest path between two points in a hydrofabric network
#' @param start_id the starting NHDPlusV2 COMID id of the shortest path. datatype: int / vector of int e.g., 61297116 or c(61297116 , 6129261) 
#' @param end_id the ending NHDPlusV2 COMID id of the shortest path. datatype: int / vector of int e.g., 61297116 or c(61297116 , 6129261) 
#' @param gpkg a local gpkg file
#' @param filename If filename is provided, data will be written using the filename
#' @param lyrs layers to extract
#' @return An sf object containing the shortest path
#' @details This function identifies the shortest path between two nodes in a hydrofabric network.
#' @author Caleb Novinger  <Caleb.D.Novinger@rtx.com>
#' @author Tadd Bindas <Tadd.N.Bindas@rtx.com>
#' @export
#' @importFrom sf st_read st_write
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom dplyr filter mutate pull %>%
#' @importFrom igraph as.igraph V E shortest_paths
#' 

get_shortest_path <- function(
    start_id,
    end_id,
    gpkg = NULL,
    filename = NULL,
    lyrs = c("divides", "flowpaths", "network", "nexus", "flowpath-attributes")
  ) {
  
    local_subset <- get_subset(
      comid=end_id, 
      gpkg=gpkg,
      lyrs=lyrs
    )
    
    network <- local_subset$network

    start_node <- network %>% 
        filter(hf_id == start_id) %>% 
        pull(id) %>%
        .[1] 
    end_node <- network %>% 
        filter(hf_id == end_id) %>% 
        pull(id) %>%
        .[1] 
    
    weights <- network$lengthkm
    weights[is.na(weights)] <- 1
    
    topo_graph <- igraph::graph_from_data_frame(
      d = data.frame(from = network$id, to = network$toid, weight = weights),
      directed = TRUE
    )
    
    all_nodes <- igraph::V(topo_graph)$name  
    if (!(start_node %in% all_nodes)) {
      stop(paste("ERROR: Start node", start_id, "not found in the graph!"))
    }
    if (!(end_node %in% all_nodes)) {
      stop(paste("ERROR: End node", end_id, "not found in the graph!"))
    }
    
    shortest_path <- igraph::shortest_paths(
      topo_graph,
      from = which(igraph::V(topo_graph)$name == start_node),
      to = which(igraph::V(topo_graph)$name == end_node),
      weights = igraph::E(topo_graph)$weight,
      mode = "out"
    )
    
    path_nodes <- igraph::V(topo_graph)[shortest_path$vpath[[1]]]$name
    
    result <- list()
    
    if ("flowpaths" %in% lyrs || length(lyrs) == 0) {
      flowpaths <- local_subset$flowpaths
      flowpaths$id <- as.character(flowpaths$id)
      shortest_flowpaths <- flowpaths[flowpaths$id %in% path_nodes, ]
      result$flowpaths <- shortest_flowpaths
    }
    
    if ("divides" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        divides <- local_subset$divides
        divides$id <- as.character(divides$id)
        shortest_divides <- divides[divides$id %in% shortest_flowpaths$id, ]
        result$divides <- shortest_divides
      }
    }
    
    if ("nexus" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        nexus <- local_subset$nexus
        nexus$id <- as.character(nexus$id)
        shortest_nexus <- nexus[nexus$id %in% shortest_flowpaths$toid, ]
        result$nexus <- shortest_nexus
      }
    }
    
    if ("network" %in% lyrs || length(lyrs) == 0) {
      network <- local_subset$network
      network$id <- as.character(network$id)
      network$toid <- as.character(network$toid)
      shortest_network <- network[network$id %in% path_nodes | 
                                    network$toid %in% path_nodes, ]
      result$network <- shortest_network
    }
    
    if ("flowpath-attributes" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        flowpath_attributes <- local_subset$`flowpath-attributes`
        flowpath_attributes$id <- as.character(flowpath_attributes$id)
        shortest_flowpath_attributes <- flowpath_attributes[flowpath_attributes$id %in% shortest_flowpaths$id, ]
        result$`flowpath-attributes` <- shortest_flowpath_attributes
      }
    }
    
    if ("flowpath-attributes-ml" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        flowpath_attributes_ml <- local_subset$`flowpath-attributes-ml`
        flowpath_attributes_ml$id <- as.character(flowpath_attributes_ml$id)
        shortest_flowpath_attributes_ml <- flowpath_attributes_ml[flowpath_attributes_ml$id %in% shortest_flowpaths$id, ]
        result$`flowpath-attributes-ml` <- shortest_flowpath_attributes_ml
      }
    }
    
    if ("hydrolocations" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        hydrolocations <- local_subset$hydrolocations
        hydrolocations$id <- as.character(hydrolocations$id)
        shortest_hydrolocations <- hydrolocations[hydrolocations$id %in% shortest_flowpaths$id, ]
        result$hydrolocations <- shortest_hydrolocations
      }
    }
    
    if ("lakes" %in% lyrs || length(lyrs) == 0) {
      if (exists("network", inherits = FALSE)) {
        lakes <- local_subset$lakes
        shortest_lakes <- lakes[lakes$hf_id %in% network$hf_id, ]
        result$lakes <- shortest_lakes
      }
    }
    
    if ("pois" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_flowpaths", inherits = FALSE)) {
        pois <- local_subset$pois
        pois$id <- as.character(pois$id)
        shortest_pois <- pois[pois$id %in% shortest_flowpaths$id, ]
        result$pois <- shortest_pois
      }
    }
    
    if ("divide-attributes" %in% lyrs || length(lyrs) == 0) {
      if (exists("shortest_divides", inherits = FALSE)) {
        divide_attributes <- local_subset$`divide-attributes`
        divide_attributes$divide_id <- as.character(divide_attributes$divide_id)
        shortest_divide_attributes <- divide_attributes[divide_attributes$divide_id %in% shortest_divides$divide_id, ]
        result$`divide-attributes` <- shortest_divide_attributes
      }
    }
    
    if (!is.null(filename)) {
      gpkg_dir <- dirname(gpkg)
      output_filename <- paste0(filename, ".gpkg")
      output_gpkg_path <- file.path(gpkg_dir, output_filename)
      
      for (layer_name in names(result)) {
        sf::st_write(result[[layer_name]], output_gpkg_path, layer = layer_name, append = FALSE)
      }
    }
    
    return(result)
}