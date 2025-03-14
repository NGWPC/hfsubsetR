% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{query_set_id}
\alias{query_set_id}
\title{Set the identifier of a query}
\usage{
query_set_id(
  query,
  identifier,
  type = c("id", "comid", "hl_uri", "poi_id", "nldi_feature", "xy", "pt")
)
}
\arguments{
\item{query}{A `hfsubset_query` object}

\item{identifier}{Identifier value}

\item{type}{Identifier type}
}
\value{
`query` with the identifier included
}
\description{
Set the identifier of a query
}
