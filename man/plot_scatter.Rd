% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_scatter.R
\name{plot_scatter}
\alias{plot_scatter}
\title{plot_scatter}
\usage{
plot_scatter(
  expr,
  x,
  y,
  dge_genes = NULL,
  transform_method = "vst",
  filter_low_values = 1,
  size = 1.5,
  alpha = 0.7,
  title = "Pairwsie scatter point",
  xlim = NULL,
  ylim = NULL,
  ns_resampling = 1000,
  color = "grey50",
  dge_genes_color = "red2"
)
}
\arguments{
\item{expr}{data.frame of expressing values,
rownames are gene names and colnames are sample/group names, all columns are expressing values,
for better transform estimation, expr should be include all sample instead of only inculde column x and y, if possible}

\item{x}{x-axis name, corelated to one column name of expr}

\item{y}{y-axis name, corelated to one column name of expr}

\item{dge_genes}{character vector of differential changed genes, corelated to rownames of expr}

\item{transform_method}{matrix transform method, one of vst, rlog or none}

\item{filter_low_values}{numeric, filter rowMeans low value rows}

\item{size}{point size}

\item{alpha}{point transparency}

\item{title}{plot title}

\item{xlim}{two element numeric vector, restrict x axis, default: NULL}

\item{ylim}{two element numeric vector, restrict y axis, default: NULL}

\item{ns_resampling}{numeric, downsampling points}

\item{color}{character vector, map point color}

\item{dge_genes_color}{character vector, map dge_genes point color}
}
\value{
ggplot2 object
}
\description{
plot scatter plot for pairwise comparison
}
