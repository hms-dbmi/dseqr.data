% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl_drug_es.R
\name{dl_drug_es}
\alias{dl_drug_es}
\title{Download drug effect size data.}
\usage{
dl_drug_es(
  files = c("cmap_es_ind.qs", "l1000_drugs_es.qs", "l1000_genes_es.qs"),
  check = FALSE
)
}
\arguments{
\item{files}{Character vector of file names to download.}

\item{check}{Check that existing drug effect size data is loadable? Default
is FALSE.}
}
\value{
Downloads drug effect size data into package folder.
}
\description{
Download drug effect size data.
}
\examples{

dl_drug_es('example.qs')

}
