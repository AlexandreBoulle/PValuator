####################
# PValuator        #
# Alexandre BOULLE #
####################

#' @importFrom stringr str_glue str_split

cumulsum <- function(x, pv){
  cs <- sum(pv < x)
  return(cs)
}


add.pvalues <- function(list.pv, l){
  n <- 0
  inc <- 0
  l.new <- length(list.pv)
  diff <- l.new - l
  pv.add <- c()
  for (i in (l + 1):l.new){
    if (n == 2){
      n <- 0
      inc <-  inc - 1
    }
    pv.inf <- list.pv[i] - (abs(list.pv[i] - list.pv[i - (diff / 2 + inc)]) / 2)
    pv.sup <- list.pv[i] + (abs(list.pv[i] - list.pv[i - (diff / 2 + inc)]) / 2)
    pv.add <- c(pv.add, pv.inf, pv.sup)
    inc <- inc + 1
    n <- n + 1
  }
  # Add new p-values
  list.pv <- c(list.pv, pv.add)
  # Update the length of p-value list
  l <- l.new
  return(list(list.pv, l))
}

#' @title Dataframe preparation
#' @description Prepare a dataframe usable by the PValuator
#' @param data a table from which we can build a dataframe usable by the PValuator
#' @param column the column name of interest (p-value column : a character)
#' @param nbr.pv the minimum number of p-values
#' @export
prepare_data <- function(data, column, nbr.pv = 20){

  # Remove NA
  data <- data[!is.na(data[, column]), ]
  # Obtain characters usable by R
  data.pvalues <- as.numeric(gsub(',', '.', data[, column]))
  # data.pvalues <- gsub(',', '.', data[, column])
  ind.min.pvalue <- which.min(data.pvalues)
  data.pvalues <- format(data.pvalues, scientific = TRUE)

  # Recover the smallest exponent for a p-value and build a vector containing p-values of interest
  # min.pvalue.chr <- as.character(min(data.pvalues))
  min.pvalue.chr <- data.pvalues[ind.min.pvalue]
  # min.pvalue.chr <- as.character(min(as.numeric(data.pvalues)))
  exp <- str_split(min.pvalue.chr, "e-")
  nbr <- as.integer(exp[[1]][2])
  # list.exp <- seq(0, (nbr-1), 1)
  list.exp <- seq(0, nbr, 1)

  # Obtain a vector of p-values
  pvalues <- c()
  add.pvalue.inf <- c(5)
  add.pvalue.sup <- c(5)
  k <- 1
  # While we haven't quite values to show a distribution we continue this procedure
  while (length(pvalues) < nbr.pv){
    pvalues <- c()
    # For each exponent
    for (i in 1:length(list.exp)){
      pvalues <- c(pvalues, str_glue("1e-{list.exp[i]}"))
      # Add values for an exponent
      for (j in 1:length(add.pvalue.inf)){
        pvalues <- c(pvalues, str_glue("{add.pvalue.inf[j]}e-{list.exp[i]}"))
        pvalues <- c(pvalues, str_glue("{add.pvalue.sup[j]}e-{list.exp[i]}"))
      }
    }
    # Keep unique values and values <= 1
    pvalues <- unique(sort(as.numeric(pvalues)))
    pvalues <- pvalues[pvalues <= 1]
    pvalues <- pvalues[pvalues >= as.numeric(min.pvalue.chr)]
    if (k == 1){
      # Add other values if we can't show a distribution
      add.pvalue.inf <- c(add.pvalue.inf, add.pvalue.inf[length(add.pvalue.inf)] * 0.5)
      add.pvalue.sup <- c(add.pvalue.sup, add.pvalue.sup[length(add.pvalue.sup)] * 1.5)
      len.inf <- length(add.pvalue.inf)
      len.sup <- length(add.pvalue.sup)
      k <- k + 1
    }
    if (k == 2){
      add.pvalue.inf <- c(add.pvalue.inf, add.pvalue.inf[length(add.pvalue.inf)] * 0.5, add.pvalue.inf[length(add.pvalue.inf)] * 1.5)
      add.pvalue.sup <- c(add.pvalue.sup, add.pvalue.sup[length(add.pvalue.sup)] * 0.5, add.pvalue.sup[length(add.pvalue.sup)] * 1.5)
      k <- k + 1
    }
    if (k > 2){
      # Results
      results.inf <- add.pvalues(add.pvalue.inf, len.inf)
      results.sup <- add.pvalues(add.pvalue.sup, len.sup)
      # Recover new p-value lists
      add.pvalue.inf <- results.inf[[1]]
      add.pvalue.sup <- results.sup[[1]]
      # Recover previous size of p-value lists
      len.inf <- results.inf[[2]]
      len.sup <- results.sup[[2]]
      k <- k + 1
    }
  }

  # Build a dataframe usable for the PValuator
  df <- data.frame(PValue = as.numeric(pvalues))
  df$Number <- apply(df, 1, function(x) cumulsum(x, pv = as.numeric(data.pvalues)))

  return(df)
}
