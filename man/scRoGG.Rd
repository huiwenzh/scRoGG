% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRoGG.R
\name{scRoGG}
\alias{scRoGG}
\title{Generate proportional distribution based on essential genes for identifying robust gene-gene correlation in scRNA-seq data}
\usage{
scRoGG(
  dat,
  normalised = TRUE,
  filter = 0.1,
  ES_number = 1000,
  org = "hsa",
  bino = 0.05,
  details = T
)
}
\arguments{
\item{dat}{A cell*gene matrix that can be raw or normliased without logarithm.}

\item{normalised}{By default, performs normalisation via scattle package.     If changes to FALSE, the input data will be used in the following steps.}

\item{filter}{By default, genes that expressed in less than 0.1 of total cells will be removed. You may want to increase the filtering threshold with more in-depth sequencing platforms.}

\item{ES_number}{Number of the ES genes used in estimating proportional correlation distribution.}

\item{org}{Selects certain ES genes under the organism. Currently only hsa and mmu are available.}

\item{bino}{Bivariate distribution probability.}

\item{details}{By default, transformed data will be returned and can be used to visualise  correlation.}
}
\value{
A list with three components:      1. bi_zero: joint zero for each gene pair.     2. prop_distribution: correlation distribution for each gene pair.    3.  transformed_data: transformed data using the mean of essential genes.
}
\description{
Generate proportional distribution based on essential genes for identifying robust gene-gene correlation in scRNA-seq data
}
\author{
Huiwen Zheng
}
