#' diann_matrix
#'
#' Generate matrix with quantities (zero will be replaced by NA values) and eventually peptides count.
#'
#' @param x data, output from diann_load
#' @param id.header Id column name (protein, gene, ...)
#' @param quantity.header Quantity column name
#' @param proteotypic.only logical; Only proteotypic peptides and the respective proteins should be considered ?
#' @param q Precursor q-value threshold
#' @param protein.q Uniquely identified protein q-value threshold
#' @param pg.q Protein group q-value threshold
#' @param gg.q Gene group q-value threshold
#' @param get_pep logical; get peptide count ?
#' @param only_pepall logical; should only keep peptide counts all or also peptide counts for each fractions ?
#'
#' @return A dataframe containing the quantities from the id you selected
#'
#' @export

diann_matrix <- function (x, id.header = "Precursor.Id", quantity.header = "Precursor.Normalised",
                          proteotypic.only = F, q = 0.01, protein.q = 1, pg.q = 1,
                          gg.q = 1, get_pep = FALSE, only_pepall = FALSE){
  df <- data.table::as.data.table(x)
  if (proteotypic.only)
    df <- df[which(df[["Proteotypic"]] != 0), ]

  dft <- df[which(df[[id.header]] != "" & df[[quantity.header]] >
                    0 & df[["Q.Value"]] <= q & df[["Protein.Q.Value"]] <=
                    protein.q & df[["PG.Q.Value"]] <= pg.q & df[["GG.Q.Value"]] <=
                    gg.q),]
  df <- unique(dft[, c("File.Name", id.header, quantity.header), with = FALSE])

  is_duplicated = any(duplicated(paste0(df[["File.Name"]],
                                        ":", df[[id.header]])))
  if (is_duplicated) {
    warning("Multiple quantities per id: the maximum of these will be calculated")
    out <- pivot_aggregate(df, "File.Name", id.header, quantity.header)
  }
  else {
    out <- pivot(df, "File.Name", id.header, quantity.header)
  }

  if(get_pep){
    x <- dft[,c("File.Name", id.header, "Genes.MaxLFQ.Unique", "Precursor.Id"), with = FALSE]
    pep <- as.data.frame(matrix(0, nrow = nrow(out), ncol = 2 + length(unique(x$File.Name))))
    colnames(pep) <- c(id.header, "peptides_counts_all", paste0("pep_count_", unique(x$File.Name)))
    for (i in unique(x$File.Name)){
      frac <- which(x$File.Name == i)
      a <- x[frac,]
      r = 1
      for (k in unique(df[[id.header]])){
        gene <- which(a[[id.header]] == k)
        b <- a[gene,]
        np <- length(unique(b[["Precursor.Id"]]))
        pep[r, c(id.header, paste0("pep_count_", i))] <- c(k, np)
        r = r + 1
      }
    }
    g <- pep[[id.header]]
    pep[[id.header]] <- NULL
    pep <- as.data.frame(apply(pep, 2, as.numeric))
    pep$peptides_counts_all <- apply(pep, 1, max)
    rownames(pep) <- g

    if(only_pepall){
      pep <- pep[,1, drop = FALSE]
    }

    out <- as.data.frame(out)
    out <- merge(out, pep, by="row.names")
    rownames(out) <- out$Row.names
    out$Row.names <- NULL
  }

  return(out)
}


### interns functions from diann-rpackage from vdemichev

pivot <- function (df, sample.header, id.header, quantity.header){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- as.data.frame(data.table::dcast.data.table(x, as.formula(paste0(id.header,
                                                                         "~", sample.header)),
                                                    value.var = "value"))
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  if(ncol(piv) == 1){
    piv$V <- 1:nrow(piv)
    piv <- piv[order(rownames(piv)), ]
    piv$V <- NULL
  }
  else{
    piv <- piv[order(rownames(piv)), ]
  }
  as.matrix(piv)
}

pivot_aggregate <- function (df, sample.header, id.header, quantity.header){
  x <- data.table::melt.data.table(df, id.vars = c(sample.header, id.header),
                                   measure.vars = c(quantity.header))
  x$value[which(x$value == 0)] <- NA
  piv <- as.data.frame(data.table::dcast.data.table(x, as.formula(paste0(id.header,
                                                                         "~", sample.header)),
                                                    value.var = "value",
                                                    fun.aggregate = function(x) max(x,
                                                                                    na.rm = TRUE)))
  rownames(piv) <- piv[[1]]
  piv[[1]] <- NULL
  if(ncol(piv) == 1){
    piv$V <- 1:nrow(piv)
    piv <- piv[order(rownames(piv)), ]
    piv$V <- NULL
  }
  else{
    piv <- piv[order(rownames(piv)), ]
  }
  piv = as.matrix(piv)
  piv[is.infinite(piv)] <- NA
  piv
}


