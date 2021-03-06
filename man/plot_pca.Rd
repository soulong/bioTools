% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pca.R
\name{plot_pca}
\alias{plot_pca}
\title{plot pca}
\usage{
plot_pca(
  data,
  metadata = NULL,
  color_by = NULL,
  shape_by = NULL,
  center = TRUE,
  scale = FALSE,
  return_data = FALSE
)
}
\arguments{
\item{data}{A data-matrix or data-frame containing numerical data only,
variables are expected to be in the rows and samples in the columns by default}

\item{metadata}{data-frame containing metadata,
strictly enforced that rownames(metadata) == colnames(mat). DEFAULT: NULL}

\item{color_by}{color point by, one of column in metadata}

\item{shape_by}{shape point by, one of column in metadata}

\item{center}{center the data before performing PCA, same as prcomp() 'center' parameter}

\item{scale}{scale the data before performing PCA, same as prcomp() 'scale' parameter}
}
\value{
list of ggplot2 object, or list of pca object
}
\description{
pca plot for a matrix
}
