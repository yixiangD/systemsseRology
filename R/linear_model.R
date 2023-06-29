# roll back

# library(lmtest)
cal_volcano_simple <-
  function(data,
           colstart,
           colend,
           null_list,
           add_variable,
           like) {
    like <- 1
    X <- as.matrix(as.data.frame(lapply(data, as.numeric)))
    # X[,colstart:ncol(X)] <-scale(X[,colstart:ncol(X)], center = T, scale = T)
    X[, 1:ncol(X)] <- scale(X[, 1:ncol(X)], center = T, scale = T)
    data.z <- as.data.frame(X)
    # name of all luminex data
    feature.name <- colnames(data.z[, colstart:colend])
    # name of demographic data149-12
    ndemo <- 1
    # matrix with rows (lumninex), cols(demographic)
    x <- rep(0, length(feature.name))
    rs <- matrix(rep(x, ndemo), ncol = ndemo)
    # save rsquare
    rsquared_terms <- data.frame(y = feature.name, rs)
    colnames(rsquared_terms) <- c("y", add_variable)
    # save pvalue
    pvalue_all <- data.frame(y = feature.name, rs)
    colnames(pvalue_all) <- c("y", add_variable)
    # save coef
    beta_all <- data.frame(y = feature.name, rs)
    colnames(beta_all) <- c("y", add_variable)
    # save tvalue
    t.value_all <- data.frame(y = feature.name, rs)
    colnames(t.value_all) <- c("y", add_variable)
    kk <- 0
    nfeat <- colend - colstart + 1
    if (like == 1) {
      # not single, likelihood ratio test
      for (feat.index in colstart:colend) {
        kk <- kk + 1
        data.z.y <- data.z
        # current luminex data
        y <- feature.name[kk]
        print(y)
        # change column name of current luminex data to "y"
        colnames(data.z.y)[feat.index] <- "y"
        all_list <- c(add_variable, null_list)
        # formula for null and alternative model
        nullformula <-
          paste("y~ ", paste(null_list, collapse = " + "))
        print(paste("null: ", nullformula))
        allformula <-
          paste("y~ ", paste(all_list, collapse = " + "))
        print(paste("all: ", allformula))
        # null and alternative model
        model.null <- lm(nullformula, data = data.z.y)
        modelsum <- summary(model.null)
        model.alternative <- lm(allformula, data = data.z.y)
        modelsum_alter <- summary(model.alternative)
        # rsquare
        rsquared_terms[(kk), 2] <- modelsum_alter[["r.squared"]]
        # values for volcanoplot
        LRT <- lrtest(model.null, model.alternative)
        pvalue <- LRT$`Pr(>Chisq)`[2]
        beta <- modelsum_alter$coefficients[add_variable, 1]
        t.value <- modelsum_alter$coefficients[add_variable, 3]
        pvalue_all[(kk), 2] <- pvalue
        beta_all[(kk), 2] <- beta
        t.value_all[(kk), 2] <- t.value
      } ## luminex
    } else {
      # single
      for (feat.index in colstart:colend) {
        kk <- kk + 1
        data.z.y <- data.z
        # current luminex data
        y <- feature.name[kk]
        print(y)
        # change column name of current luminex data to "y"
        colnames(data.z.y)[feat.index] <- "y"
        # formula for null and alternative model
        nullformula <-
          paste("y~ ", paste(null_list, collapse = " + "))
        print(paste("null: ", nullformula))
        # null and alternative model
        model.null <- lm(nullformula, data = data.z.y)
        modelsum <- summary(model.null)
        # # logistic model
        # model.null <- glm("pf ~ y", data = data.z.y, family= "binomial")
        # modelsum<-summary(model.null)
        # rsquare
        rsquared_terms[(kk), 2] <- modelsum[["r.squared"]]
        # values for volcanoplot
        beta <- modelsum$coefficients[add_variable, 1]
        t.value <- modelsum$coefficients[add_variable, 3]
        pvalue <- modelsum$coefficients[add_variable, 4]
        # beta <- modelsum$coefficients["y", 1]
        # t.value <- modelsum$coefficients["y", 3]
        # pvalue <- modelsum$coefficients["y", 4]
        #
        pvalue_all[(kk), 2] <- pvalue
        beta_all[(kk), 2] <- beta
        t.value_all[(kk), 2] <- t.value
      } ## luminex
    }
    # for all demographic data
    kk <- 1
    r2 <- rsquared_terms[, kk + 1]
    t.value <- t.value_all[, kk + 1]
    beta <- beta_all[, kk + 1]
    pvalue <- pvalue_all[, kk + 1]
    labels <- feature.name
    tb <- rep("High", length(t.value))
    tb[which(t.value < 0)] <- "Low"
    pvalue <- p.adjust(pvalue, method = "BH", n = length(pvalue))
    lmm.df <-
      data.frame(
        feature = feature.name,
        tvalue = t.value,
        coefficient = beta,
        lab = labels,
        assay = gsub("_.*", "", labels),
        group = tb,
        r2 = r2,
        pvalue = pvalue
      )
    lmm.save <-
      data.frame(
        group = tb,
        feature = labels,
        tvalue = t.value,
        pvalue = pvalue
      )
    lmm.df5 <- lmm.df[which(lmm.df$r2 > 0.3), ]
    if (like == 1) {
      lmm <- list(lmm.df, nullformula, allformula)
    } else {
      lmm <- list(lmm.df, nullformula)
    }
    return(lmm)
  }


plot_volcano <-
  function(data,
           filename,
           title,
           foldername,
           pthreshold,
           col.arms) {
    # library(ggrepel)
    setwd(foldername)

    g.arms <- c("High", "Low")
    # col.arms <- c("#ffb300", "#b300ff")
    #
    p <-
      ggplot(data, aes(coefficient, -log10(pvalue), label = lab)) +
      geom_hline(
        yintercept = -log10(pthreshold),
        colour = "red",
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = 0,
        colour = "slategray",
        linetype = "dashed"
      ) +
      scale_fill_manual(breaks = g.arms, values = col.arms) +
      geom_text_repel(
        size = 5.5,
        aes(colour = assay),
        show.legend = FALSE,
        max.overlaps = 30
      ) +
      geom_point(size = 5, aes(fill = group), shape = 21) +
      theme_classic() +
      xlab("Coefficient") +
      ylab("-log10(p-value)") +
      ggtitle(title) +
      theme(
        text = element_text(size = 16.2),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 24)
      ) +
      theme(legend.title = element_blank(), legend.position = "none")
    print(p)
    ggsave(
      filename = paste0(filename, ".pdf"),
      p,
      width = 10,
      height = 10,
      dpi = 300,
      device = "pdf"
    )
    return(p)
  }
