#' Find an origin from indexed IDs
#' @param network A `dplyr`-compatible object.
#' @param id A queryable identifier of type `type`.
#' @param type An index type describing `id`.
#' @returns A network origin. If a single origin is not found,
#'          then an exception is raised.
#' @export
find_origin <- function(
  network,
  id,
  src,
  type = c("id", "comid", "hl_uri", "poi_id", "nldi_feature", "xy")
) {

  hydroseq <- NULL
  type <- match.arg(type)
  query <- structure(id, class = type)

  origin <- try(find_origin_query(query, network, src))

  if (inherits(origin, "try-error")) {
    stop(origin, call. = FALSE)
  }
  
  origin <-
    origin |>
    dplyr::select(any_of(c('id', 'toid', 'vpuid', 'topo', 'hydroseq'))) |>
    dplyr::distinct() |>
    dplyr::collect() |> 
    dplyr::slice_min(hydroseq, with_ties = TRUE)

  if (nrow(origin) == 0) {
    stop("No origin found")
  } else if (nrow(origin) > 1) {
    stop("Multiple origins found: ", dput(origin$id))
  } else {
    return(origin)
  }
}

#' S3 dispatch on query identifier type
#' @param id A queryable identifier, see `find_origin`.
#' @param network A `dplyr`-compatible object.
#' @param src gpkg layer
#' @returns `network` after applying a [dplyr::filter] expression.
#' @keywords internal
find_origin_query <- function(id, network, src) {
  UseMethod("find_origin_query")
}

#' @method find_origin_query default
#' @keywords internal
find_origin_query.default <- function(id, network, src) {
  stop(paste(
    "identifier of class",
    paste0("`", class(id), "`", collapse = "/"),
    "not supported"
  ))
}

#' @method find_origin_query id
#' @keywords internal
find_origin_query.id <- function(id, network, src) {
  id <- unclass(id)
  dplyr::filter(network, id == !!id)
}

#' @method find_origin_query comid
#' @keywords internal
find_origin_query.comid <- function(comid, network, src) {
  hf_id <- NULL
  comid <- unclass(comid)
  dplyr::filter(network, hf_id == !!comid)
}

#' @method find_origin_query hl_uri
#' @keywords internal
find_origin_query.hl_uri <- function(hl_uri, network, src) {
  hl_uri <- unclass(hl_uri)
  dplyr::filter(network, hl_uri == !!hl_uri)
}

#' @method find_origin_query poi_id
#' @keywords internal
find_origin_query.poi_id <- function(poi_id, network, src) {
  poi_id <- unclass(poi_id)
  dplyr::filter(network, poi_id == !!poi_id)
}

#' @method find_origin_query nldi_feature
#' @keywords internal
find_origin_query.nldi_feature <- function(nldi_feature, network, src) {
  .Class <- "comid"

  nldi_feature <- structure(
    nhdplusTools::discover_nhdplus_id(nldi_feature = nldi_feature),
    class = "comid"
  )

  NextMethod()
}

#' @method find_origin_query xy
#' @keywords internal
find_origin_query.xy <- function(xy, network, src) {
  .Class <- "id"

  if(grepl("https", src)){
    src = paste0("/vsicurl/", src)
  } else if(grepl("s3", src)){
    src = paste0("/vsis3/", src)
  } else {
    src = src
  }
  
  crs <- as_ogr(src, "divides") |> 
    head(1) |> 
    sf::st_as_sf() |> 
    st_crs()
  
  bb <- sf::st_point(unclass(xy)) |> 
    sf::st_sfc(crs = 4326) |> 
    sf::st_as_sf() |> 
    sf::st_transform(crs$wkt) |> 
    sf::st_geometry() |> 
    sf::st_as_text()

  xy <- structure(
    sf::read_sf(src, "divides", wkt_filter = bb) |> 
      dplyr::pull(id),
    class = "id"
  )

  NextMethod()
}
