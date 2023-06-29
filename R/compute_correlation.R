#' Compute correlation for given columns A and columns B
#'
#' @param data n_samples x n_features matrix
#' @param ab.labels features group A
#' @param tcell.labels features group B
#'
#' @return df data.frame with correlation and p value
#' @export
#'
get_cor_col <- function(data, ab.labels, tcell.labels) {
  res <-
    Hmisc::rcorr(as.matrix(data[, c(ab.labels, tcell.labels)]), type = "spearman")
  resr <- res$r
  resr <- as.data.frame(resr[ab.labels, tcell.labels])
  resr.mat <- as.matrix(resr)
  resr$antibody <- row.names(resr)

  resp <- res$P
  resp <- as.data.frame(resp[ab.labels, tcell.labels])
  resp$antibody <- row.names(resp)

  df.r <- resr %>% tidyr::gather(tcell, value, -c(antibody))
  df.p <- resp %>% tidyr::gather(tcell, p, -c(antibody))
  df.p$q <- p.adjust(df.p$p, method = "BH")

  identical(df.r$antibody, df.p$antibody)
  identical(df.r$tcell, df.p$tcell)

  df <- cbind(df.r, df.p[, c("p", "q")])
  colnames(df)[1:3] <- c("from", "to", "cor")
  return(df)
}
