% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tx2gene.R
\name{build_ensdb}
\alias{build_ensdb}
\title{Build and install ensembldb annotation package.}
\usage{
build_ensdb(species = "Homo sapiens", release = "94")
}
\arguments{
\item{species}{Character vector indicating species. Genus and species should
be space separated, not underscore. Default is \code{'Homo sapiens'}.}

\item{release}{EnsemblDB release. Should be same as used in
\code{\link[rkal]{build_kallisto_index}}.}
}
\value{
Called for side effects.
}
\description{
Build and install ensembldb annotation package.
}
\examples{

# build ensembldb annotation package for human
# build_ensdb()
}
\keyword{internal}
