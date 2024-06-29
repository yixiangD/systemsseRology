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
    return(list(feat = sel_features, df_score = NULL))
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
    ggsave(paste(saved.dir, paste(prefix, "score.pdf", sep = "_"), sep = "/"), plt_scores + theme(legend.position = "top"), width = 4, height = 3)
    opts_plot$LV_ind <- 1
    lvs <- ropls::getScoreMN(model)
    df_scores <- data.frame(lvs)
    df_scores$y <- y
    colnames(df_scores) <- c("LV1", "LV2", "group")

    # set additional options required to color code enrichment in the bar plot of the loadings
    plt_loadings_bar <-
      my_visualize_ropls_loadings_bar(model, options = opts_plot)
    ggsave(
      paste(saved.dir, paste(prefix, "lv1.pdf", sep = "_"), sep = "/"),
      plt_loadings_bar + theme(legend.position = "top"),
      width = 4,
      height = 3
    )
    opts_plot$LV_ind <- 2
    plt_loadings_bar2 <-
      my_visualize_ropls_loadings_bar(model, options = opts_plot)
    ggsave(
      paste(saved.dir, paste(prefix, "lv2.pdf", sep = "_"), sep = "/"),
      plt_loadings_bar2 + theme(legend.position = "top"),
      width = 4,
      height = 3
    )
    return(list(feat = sel_features, df_score = df_scores))
  }
}


my_visualize_ropls_loadings_bar <- function(model, options = list()) {
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
  # according to enrichent in classes
  if (options$mark_enrichment) {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV, fill = mark)) +
      ggplot2::scale_fill_manual(values = options$colors[[y_name]])
  } else {
    plt_bar <- ggplot2::ggplot(data = df_loadings, ggplot2::aes(x = features, y = LV))
  }
  plt_bar <- plt_bar +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab(paste("LV", options$LV_ind, " loadings", sep = "")) +
    ggplot2::labs(fill = "enriched in") +
    ggplot2::scale_x_discrete(labels = df_loadings$labels) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(color = "black")) # ,
  # axis.text.y = element_text(colour = as.character(feature_annot$useColor[match(dfBar$features[order(dfBar$vipScores)],
  # rownames(feature_annot))])))
}
