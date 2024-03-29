% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backup.R
\name{getProbe}
\alias{getProbe}
\title{Regional Plot Association Results of Pig}
\usage{
getProbe(geneName, GPL, expr)
}
\arguments{
\item{geneName}{a character of gene SYMBOL name}

\item{GPL}{a character vector of GPL number}

\item{expr}{a GSE expression dataframe where probes would to be selected based on their probes max expression}
}
\value{
a character of probe
}
\description{
This is a plotting tool to simultaneously display the association p-value,
             the LD, recombination rate and annotated gene in region of interesting in pig.
             The interpretation and visualization tool of association results will
             contribute to further identification of causal mutations.
}
\details{
nothing
}
