% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMean.R
\name{getMean}
\alias{getMean}
\title{Run Linear Models}
\usage{
getMean(x)
}
\arguments{
\item{x}{A phyloseq object which combined OTU count, taxonomy and metadata}
}
\value{
A matrix of the infile
}
\description{
This function takes in phyloseq object, association matix, p value matrix
and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
define OTU-Phontype assocaition.
}
