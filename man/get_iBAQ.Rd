% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_iBAQ.R
\name{get_iBAQ}
\alias{get_iBAQ}
\title{get_iBAQ}
\usage{
get_iBAQ(
  dataset,
  proteinDB = list(),
  id_name = "Protein.Group",
  ecol = 5:7,
  peptideLength = c(5, 36),
  nbMiscleavages = 0,
  proteaseRegExp = getProtease("trypsin"),
  log2_transformed = TRUE,
  keep_original = TRUE
)
}
\arguments{
\item{dataset}{A data frame that you want to get the iBAQ quantities from.}

\item{proteinDB}{A list that contains for each protein the theoretical peptide sequence.
It doesn't have to be the same length as the number of proteins from the dataset;
if it's the case, value from proteins that are not in this list will be set to NA.}

\item{id_name}{The name of the column that contains the protein IDs.}

\item{ecol}{The indexes of the intensity columns.}

\item{peptideLength}{A vector that contains the minimum and maximum peptide length.}

\item{nbMiscleavages}{The number of miss cleavages.}

\item{proteaseRegExp}{The enzyme you used to digest proteins.}

\item{log2_transformed}{Logical to tell if you want to log2 transformed the data.}

\item{keep_original}{Logical to tell if you want to keep the original intensity values.}
}
\value{
The dataset with the iBAQ quantities and number of theoretical peptides.
}
\description{
Function to calculate iBAQ quantities from a dataset.
}
\details{
This function use three function from \pkg{SafeQuant} package from github eahrne/SafeQuant. Those three function
are simple and small, so to facilitate dependencies and lower the size of this package, I chose to copy-paste those three.

Erik Ahrne (2021). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.4.
}
\examples{
library(DIAgui)
data <- small_report \%>\% dplyr::filter(Q.Value <= 0.01 & PG.Q.Value <= 0.01)
data <- diann_maxlfq(data,
                  group.header="Protein.Group",
                  id.header = "Precursor.Id",
                  quantity.header = "Precursor.Normalised"
                  )
nc <- ncol(data)
data$Protein.Group <- rownames(data)
rownames(data) <- 1:nrow(data)
data <- data[order(data$Protein.Group),]
data <- data[,c((nc+1):ncol(data), 1:nc)]
sequence_list <- getallseq(pr_id = data$Protein.Group)
iBAQ <- get_iBAQ(data, sequence_list, ecol = 2:4)

}
