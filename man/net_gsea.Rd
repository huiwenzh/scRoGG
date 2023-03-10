% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net_gsea.R
\name{net_gsea}
\alias{net_gsea}
\title{Perform network-based geneset enrichment analysis}
\usage{
net_gsea(G, species='Homo sapiens',category='H', minSize=5)
}
\arguments{
\item{species}{Species information to msigdbr package.}

\item{category}{Category information to msigdbr package.}

\item{subcategory}{Subcategory information to msigdbr package.}

\item{minSize}{Minimum size in each pathway.}

\item{maxSize}{Maximum size in each pathway.}

\item{BPPARAM}{Parellalisation}
}
\description{
Perform network-based geneset enrichment analysis
}
\author{
Huiwen Zheng
}
