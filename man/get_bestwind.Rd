% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bestwind.R
\name{get_bestwind}
\alias{get_bestwind}
\title{get_bestwind}
\usage{
get_bestwind(
  data,
  n_window = 25,
  per_frac = FALSE,
  overlap = 0,
  window_size = NULL
)
}
\arguments{
\item{data}{The report-lib data from DIA-nn which contains the \emph{PrecursorMz} column.
It can be a path to the file or the dataframe corresponding to the file
(need to have \emph{FileName} and \emph{PrecursorMz} columns in that case).}

\item{n_window}{The number of windows you want to have.}

\item{per_frac}{If \code{FALSE}, will select the best windows from all fraction without differentiation.
Else, it will select the best window for each fraction and then do the mean of those.}

\item{overlap}{A fixed overlap between windows. Default is 0.}

\item{window_size}{A fix m/z window size. Default is NULL but if numeric, it will compute n
windows of same size.}
}
\value{
A list containing a data frame with \emph{FileName}, \emph{PrecursorMz} and a \emph{bins} column;
        four plots corresponding to bar plots / histogramm from the distribution of the \emph{PrecursorMz}
        with and without the window selection.
}
\description{
Get the best m/z windows for DIA according a number of window and the report-lib from DIA-nn.
}
\details{
In the end, the distribution of the \emph{PrecursorMz} per window will be uniform.
}
