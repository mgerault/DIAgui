% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcanoDIA.R
\name{volcanoDIA}
\alias{volcanoDIA}
\title{volcanoDIA}
\usage{
volcanoDIA(
  data,
  control,
  treated,
  id = NULL,
  transformation = c("none", "log2"),
  tit = "",
  FDR = 0.01,
  FC_cut = 2.5,
  curvature = 0.1,
  save_file = FALSE
)
}
\arguments{
\item{data}{Processed data from DIAnn (from iq processing or diann_matrix for example)}

\item{control}{The names of the controls in your data (needs to match column names in your data)}

\item{treated}{The names of the treated in your data (needs to match column names in your data)}

\item{id}{A character indicating the name of the coulumns that conatins the ids.
If NULL, will check for the column named 'id' or the rownames.}

\item{transformation}{Which transformation do you want to apply (log2 or none)}

\item{tit}{Title of your plot}

\item{FDR}{The FDR for the volcano plot}

\item{FC_cut}{The fold-change cutoff}

\item{curvature}{The curvature used for the curve on the volcano plot}

\item{save_file}{Save results file}
}
\value{
A volcano plot
}
\description{
Plot volcano plot from your DIA data
}
