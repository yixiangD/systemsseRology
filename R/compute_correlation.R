#' Compute correlation for given columns A and columns B
#'
#' @param data: n_samples x n_features matrix
#' @param from.labels: features group A, vector of strings
#' @param to.labels: features group B, vector of strings
#'
#' @return df data.frame with correlation and p value
#' @export
#'
get_cor_col <- function(data, from.labels, to.labels) {
  res <-
    Hmisc::rcorr(as.matrix(data[, c(from.labels, to.labels)]), type = "spearman")
  resr <- res$r
  resr <- as.data.frame(resr[from.labels, to.labels])

  resp <- res$P
  resp <- as.data.frame(resp[from.labels, to.labels])
  # create a new column with name "from"
  resr$from <- row.names(resr)
  resp$from <- row.names(resp)
  rownames(resr) <- NULL
  rownames(resp) <- NULL
  df.r <- resr %>% tidyr::pivot_longer(cols = !from, names_to = "to", values_to = "cor")
  df.p <- resp %>% tidyr::pivot_longer(cols = !from, names_to = "to", values_to = "p")
  df.p$q <- p.adjust(df.p$p, method = "BH")
  # test exact equality
  # identical(df.r$from, df.p$from)
  # identical(df.r$to, df.p$to)
  df <- merge(df.r, df.p, by = c("from", "to"))
  return(df)
}
# if (!keep.tri %in% c("U", "L", "A")) {
#  stop("please provide keep.tri within U, L or A")
# }else if (keep.tri == "U") {
#    resr[upper.tri(resr)] <- NA
# } else if (keep.tri == "L") {
#    resr[upper.tri(resr)] <- NA
# }
