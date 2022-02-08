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
#' @param Top3 logical; get Top3 absolute quantification
#'
#' @return A dataframe containing the quantities from the id you selected
#'
#' @export

diann_matrix <- function (x, id.header = "Precursor.Id", quantity.header = "Precursor.Normalised",
                          proteotypic.only = FALSE, q = 0.01, protein.q = 1, pg.q = 1,
                          gg.q = 1, get_pep = FALSE, only_pepall = FALSE, margin = -10, Top3 = FALSE){
  df <- data.table::as.data.table(x)
  if(proteotypic.only){
    df <- df[which(df[["Proteotypic"]] != 0), ]
  }
  dft <- df[which(df[[id.header]] != "" & df[["Q.Value"]] <= q & df[["Protein.Q.Value"]] <=
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
  if(Top3 & id.header == "Genes"){
    x <- unique(dft[which(dft[["Genes"]] != ""), c("File.Name",
                                                   "Genes", "Precursor.Id", "Precursor.Normalised"), with = FALSE])
    x[["File.Name"]] <- as.character(x[["File.Name"]])
    x[["Genes"]] <- as.character(x[["Genes"]])
    x[["Precursor.Id"]] <- as.character(x[["Precursor.Id"]])
    x[["Precursor.Normalised"]] <- as.numeric(x[["Precursor.Normalised"]])
    if (any(x[["Precursor.Normalised"]] < 0, na.rm = T))
      stop("Only non-negative quantities accepted")
    is_duplicated = any(duplicated(paste0(x[["File.Name"]],
                                          ":", x[["Genes"]], ":", x[["Precursor.Id"]])))
    if (is_duplicated)
      warning("Multiple quantities per id: the maximum of these will be calculated")
    x[["Precursor.Normalised"]][which(x[["Precursor.Normalised"]] == 0)] <- NA
    x[["Precursor.Normalised"]] <- log(x[["Precursor.Normalised"]])
    x[["Precursor.Normalised"]][which(x[["Precursor.Normalised"]] <= margin)] <- NA
    x <- x[!is.na(x[["Precursor.Normalised"]]), ]
    genes <- unique(x[["Genes"]])
    samples <- unique(x[["File.Name"]])
    top3_res <- list()
    for(i in genes){
      top3 <- x[which(x[["Genes"]] == i),]
      top3[["Genes"]] <- NULL
      top3 <- tidyr::spread(top3, File.Name, Precursor.Normalised)
      top3[["Precursor.Id"]] <- NULL
      top3_res[[i]] <- top3
    }
    top3_res <- lapply(top3_res, function(x){
      x <- apply(x, 2, function(y){
        if(sum(!is.na(y)) < 3){
          y <- NA
        }
        else{
          y <- y[order(y, decreasing = TRUE)][1:3]
          y <- mean(y)
        };
        y
      })
      x <- as.data.frame(t(x))
      n <- samples[!(samples %in% colnames(x))]
      if(!purrr::is_empty(n)){
        for(i in n){
          x[[i]] <- NA
        }
      };
      x
    })
    p = names(top3_res)
    top3_res <- Reduce(rbind, top3_res)
    rownames(top3_res) <- p
    top3_res <- top3_res[order(rownames(top3_res)),]
    colnames(top3_res) <- paste0("Top3_", colnames(top3_res))
    top3_res <- exp(top3_res)

    out <- as.data.frame(out)
    out <- merge(out, top3_res, by="row.names")
    rownames(out) <- out$Row.names
    out$Row.names <- NULL
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


