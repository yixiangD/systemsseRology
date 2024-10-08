#' Partial Least Squares Discriminant Analysis (PLS-DA)
#' This function performs Partial Least Squares Discriminant Analysis on the given dataset.
#'
#' @param X, matrix of predictive features, rows for samples and columns for measurements
#' @param y, vector of outcome variables
#' @param clr.grps, color palettes vector, with names corresponding to possible
#' outcomes
#' @param saved.dir, exporting directory/path
#' @param prefix, prefix to the file names
#' @param feature_selection_th, feature selection threshold (optional)
#' @return a list [sel_features, df_score, df_load1, df_load2] containing selected features and
#' dataframe of latent space scores & dataframe of loadings
#' @export

pls_da <- function(X, y, clr.grps, saved.dir, prefix, feature_selection_th = 0.8) {
  opts_sel <- list(n_trials = 50, threshold = feature_selection_th, return_count = FALSE)

  opts_model <- list(n_LV = 2)

  opts_plot <- list(
    loading_alpha = 1.0, # transparency for the loadings
    score_alpha = 1.0, # transparency for the scores
    LV_ind = c(1, 2), # which LVs to plot
    colors = clr.grps,
    y_name = "group"
  )

  opts_plot$mark_enrichment <- TRUE

  opts_model <- list(n_LV = 2)

  X <- as.matrix(X)
  X <- scale(X)
  y <- as.factor(y)

  sel_features <- select_repeat(X, y, selector = select_lasso, options = opts_sel)

  # The model only has one latent variable. For visualization purposes we fix it to be two dimensional.
  if (length(sel_features) == 1) {
    return(list(feat = sel_features, df_score = NULL, df_load = NULL))
  } else {
    X_sel <- X[, sel_features]
    opts_plot$X <- X_sel
    opts_plot$y <- y
    model <- train_ropls(X_sel, y, options = opts_model)
    # add ROC
    true_labels_numeric.vi <- round(as.numeric(model@suppLs[["yModelMN"]]))
    predicted_values.vn <- as.numeric(model@suppLs[["yPreMN"]])
    sac.roc <- pROC::roc(true_labels_numeric.vi, predicted_values.vn)
    pdf(paste(saved.dir, paste(prefix, "roc.pdf", sep = "_"), sep = "/"), width = 4, height = 3)
    plot(sac.roc, print.auc = TRUE)
    dev.off()

    plt_scores <-
      visualize_ropls_scores(model, y, options = opts_plot)

    plt_scores_top_legend <- plt_scores + ggplot2::theme(legend.position = "top")
    ggplot2::ggsave(paste(saved.dir, paste(prefix, "score.pdf", sep = "_"), sep = "/"), plt_scores_top_legend, width = 4, height = 3)

    lvs <- ropls::getScoreMN(model)
    df_scores <- data.frame(lvs)
    df_scores$y <- y
    colnames(df_scores) <- c("LV1", "LV2", "group")

    # set additional options required to color code enrichment in the bar plot of the loadings
    # barplot LV1
    opts_plot$LV_ind <- 1
    df_load1 <- pls_da_bar_stats(model, options = opts_plot)

    # barplot LV2
    opts_plot$LV_ind <- 2
    df_load2 <- pls_da_bar_stats(model, options = opts_plot)
    return(list(feat = sel_features, df_score = df_scores, df_load1 = df_load1, df_load2 = df_load2))
  }
}


#' My customized plotting schemes for PLS-DA barplots
#' This function performs customized plotting for PLS-DA.
#'
#' @param model, trained ropls object
#' @param options, list of plotting options
#' @return df_load, data.frame containing loading values for each latent dimension
#' @export
pls_da_bar_stats <- function(model, options = list()) {
  # ----------------- BEGIN OPTIONS ----------------- #
  if (!("LV_ind" %in% names(options))) {
    options$LV_ind <- c(1)
  }
  if ("y" %in% names(options)) {
    y <- options$y
    if (is.factor(y)) {
      n_groups <- nlevels(y)
    } else {
      n_groups <- NA
    }
  } else {
    n_groups <- NA
  }
  if (!("mark_enrichment" %in% names(options)) | is.na(n_groups)) {
    options$mark_enrichment <- FALSE
  }
  if (options$mark_enrichment & (is.na(n_groups) | !("X" %in% names(options)))) {
    stop("Enrichment only works for classification and when X and y are provided")
  }
  # color for the scores and name of the grouping
  if (!("y_name" %in% names(options))) {
    y_name <- "y"
  } else {
    y_name <- options$y_name
  }
  # color for the scores and name of the grouping
  if (!("colors" %in% names(options)) | length(grep(y_name, names(options$colors))) == 0) {
    if (is.factor(y)) {
      tmp <- rep(NA, length = nlevels(y))
      names(tmp) <- levels(y)
      for (ind in 1:nlevels(y)) {
        tmp[ind] <- RColorBrewer::brewer.pal(n = max(3, nlevels(y)), name = "Dark2")[ind]
      }
      options$colors <- list()
      options$colors[[y_name]] <- tmp
    } else {
      # For regression, a color palette needs to be provided
      options$colors$y <- list(low = "#C7E4F9", high = "#004D7F")
    }
  }
  if (ropls::getSummaryDF(model)$pre +
    ropls::getSummaryDF(model)$ort < options$LV_ind) {
    stop("required LV exceed existing LVs")
  }
  if (!("df_features" %in% names(options))) {
    options$df_features <- data.frame(
      name = rownames(model@loadingMN),
      label = rownames(model@loadingMN)
    )
  }
  # ----------------- END OPTIONS ----------------- #
  # check first whether its a orthogonal PLS or a regular PLS
  if (ropls::getSummaryDF(model)$ort > 0) {
    stop("orthogonal PLS-DA not supported yet")
    # if (options$LV_ind[1] == 1) {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model),
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[,options$LV_ind[2] - 1])
    # } else {
    #   df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[1] - 1],
    #                             LV2 = ropls::getLoadingMN(model, orthoL = TRUE)[, options$LV_ind[2] - 1])
    # }
  } else {
    df_loadings <- data.frame(
      LV = ropls::getLoadingMN(model)[, options$LV_ind[1]],
      vip_scores = ropls::getVipVn(model)
    )
    df_loadings$features <- rownames(df_loadings)
    df_loadings$labels <- options$df_features$label[match(rownames(df_loadings), options$df_features$name)]
  }

  # TODO: catch if its an orthogonal
  if (options$mark_enrichment & !is.na(n_groups)) {
    df_loadings$mark <- NA
    X <- options$X
    for (ind_feat in 1:nrow(df_loadings)) {
      tmp_mean <- rep(NA, length = nlevels(y))
      for (ind_class in 1:nlevels(y)) {
        tmp_mean[ind_class] <- mean(X[
          which(y == levels(y)[ind_class]),
          which(colnames(X) == df_loadings$features[ind_feat])
        ])
      }
      df_loadings$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
    }
    df_loadings$mark <- factor(df_loadings$mark, levels = levels(y))
  }
  df_loadings <- df_loadings[order(df_loadings$LV), ]
  df_loadings$features <- factor(df_loadings$features, levels = unique(df_loadings$features))
  # plot loadings sorted according to the VIP score and color coding it
  return(df_loadings)
}
