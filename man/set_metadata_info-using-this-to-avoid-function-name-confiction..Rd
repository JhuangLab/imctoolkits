% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadata.R
\name{set_metadata_info using this to avoid function name confiction.}
\alias{set_metadata_info using this to avoid function name confiction.}
\alias{set_metadata_info}
\title{set_metadata_info}
\usage{
set_metadata_info(sce, sample_id, stype = "tumor")
}
\arguments{
\item{sce}{object, a singcellexperiment object.}
}
\value{
a singcellexperiment object
}
\description{
set_metadata_info
}
\examples{

sce <- set_metadata_info(sce, sample_id = "p155_pdac")
sce <- set_metadata_info(sce, sample_id = "p155_pdac", stype = "tumor")

}
