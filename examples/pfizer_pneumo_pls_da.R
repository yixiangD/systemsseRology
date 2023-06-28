library(readxl)
library(ggpubr)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)
library("writexl")
library(tools)
library(systemsseRology)
library(ggplot2)
library(patchwork)

# load data----

rm(list = ls())

color.v1 <- "#c2bbf0"
color.v2 <- "#06aed5"
color.v3 <- "#4c956c"

color.group <- c("13v_23v" = color.v1, "13v_13v" = color.v2, "23v_13v" = color.v3)

color.13v <- "#cdb4db"
color.23v <- "#ffbe0b"

saved.dir <- "~/Downloads/pfizer"
if (!dir.exists(saved.dir)) {
  dir.create((saved.dir))
}
setwd("/Users/yd973/Library/CloudStorage/OneDrive-MassGeneralBrigham/projects/pfizer")

path <- "/Users/yd973/Dropbox (Partners HealthCare)/!Pneumo round 2"

fname <- paste(path, "data/serology_clean0426.xlsx", sep = "/")
data <- read_excel(fname)

data$group <- paste(data$`vacc 1`, data$`vacc 2`, sep="_")
colnames(data) <- gsub(" ", "_", colnames(data))

col.feat <- colnames(data)[grepl("AD|Fc|Ig", colnames(data))]
col.id <- colnames(data)[!grepl("AD|Fc|Ig", colnames(data))]
# preprocess-----
df.long <- reshape2::melt(data, id=col.id)
df.sep <- separate(df.long, col="variable", into=c("assay", "antigen"), sep="_")
df.sep$antigen <- ifelse(df.sep$antigen == "19A", "Pn19A", df.sep$antigen)
df.sep$variable <- paste(df.sep$assay, df.sep$antigen, sep="_")
col.ids <-
  c(
    "study",
    "patient_ID",
    "visit",
    "vacc_1",
    "vacc_2",
    "sex",
    "age_at_dose1",
    "group",
    "variable"
  )
df <- pivot_wider(df.sep[, c(col.ids, "value")], names_from = "variable")

df.3005 <- df[df$study == "3005",]
df.3010 <- df[df$study == "3010",]

# PLSDA -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
study <- args[1]
visit <- args[2]

# for testing-----
# study <- "3005"
# visit <- "VISIT 5"

stopifnot(study %in% c("3005", "3010"))
stopifnot(visit %in% c("VISIT 2", "VISIT 5"))

col.feat <- colnames(df)[grepl("Ig|Fc|AD", colnames(df))]

if (study == "3005") {
  df.local <- df.3005[df.3005$visit == visit, ]
} else {
  df.local <- df.3010[df.3010$visit == visit, ]
}
if (visit == "VISIT 2") {
  df.local$Vaccine <- df.local$vacc_1

  my_colors <- list(group = c(
    "13v" = color.13v,
    "23v" = color.23v
  ))
  col.feat <- col.feat[!grepl("6A", col.feat)]
} else {
  df.local$Vaccine <- df.local$group
  my_colors <- list(group = color.group)
}

col.sel <- col.feat[!grepl("CRM|HSA|HA", col.feat)]

X <- df.local[, col.sel]
# in case its a character data frame
X <- as.data.frame(lapply(X, as.numeric))
X <- as.matrix(X)
print(paste(dim(X)[2], " features in total for PLS-DA", sep = " "))

# check na in X
#X <- knnImputation(X)
X <- as.matrix(glmnet::makeX(as.data.frame(X), na.impute = TRUE))

log.ind <- which(grepl("ADCD|IgG|IgA|IgM|FcR", colnames(X)))
X[, log.ind] <- log10(X[, log.ind])
X <- scale(X, center = TRUE, scale = TRUE)

y <- df.local$Vaccine
y <- factor(y)

df_features <- data.frame(
  name = colnames(X),
  label = colnames(X)
)

# Feature selection-----

opts_sel <- list(threshold = 0.8, n_trials = 100, return_count = TRUE)
if (min(table(y)) <= 5) {
  X <- rbind(X, X)
  y <- c(y, y)
  y <- factor(y)
}

out <- select_repeat(X, y, selector = select_lasso, options = opts_sel)

df_count <- data.frame(
  features = names(out$feature_count),
  name = names(out$feature_count),
  selected = out$feature_count * 100 / opts_sel$n_trials,
  mark = NA
)

df_count <- df_count[which(df_count$selected > 0), ]
df_count <- df_count[order(-df_count$selected), ]
df_count$features <- df_features$label[match(df_count$features, df_features$name)]
df_count$features <- factor(df_count$features, levels = df_count$features)
# annotation where feature is enriched
for (ind_feat in 1:nrow(df_count)) {
  tmp_mean <- rep(NA, length = nlevels(y))
  for (ind_class in 1:nlevels(y)) {
    tmp_mean[ind_class] <- mean(X[
      which(y == levels(y)[ind_class]),
      which(colnames(X) == df_count$name[ind_feat])
    ])
  }
  df_count$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
}
df_count$mark <- factor(df_count$mark, levels = levels(y))

plt_bar <- ggplot(data = df_count, aes(x = features, y = selected, fill = mark)) +
  scale_fill_manual(values = my_colors$group) +
  geom_bar(stat = "identity", color = "black") +
  xlab("") +
  geom_hline(yintercept = opts_sel$threshold * 100) +
  ylab("selected (%)") +
  labs(fill = "enriched in") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14), axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 10, angle = 80, hjust = 1, vjust = 1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
# print(plt_bar)

sel_features <- out$sel_features

X_sel <- X[, sel_features]
opts_model <- list(n_LV = 2)
model <- train_ropls(X_sel, y, options = opts_model)

opts_plot <- list(
  LV_ind = c(1, 2), # which LVs to plot
  colors = my_colors,
  size = 2.5,
  stroke = 0.5,
  alpha = 1.0,
  y_name = "group",
  level = 0.95
)

p.scores <- visualize_ropls_scores(model, y, options = opts_plot) +
  theme(legend.title = element_blank(), legend.position = "top")

# print(p.scores)

pdf(paste(saved.dir, paste(study, visit, "pls-da_scores.pdf", sep = "_"), sep = "/"), width = 4, height = 3)
print(p.scores)
dev.off()

loadsm <- as.data.frame(model@loadingMN)
vips <- as.data.frame(model@vipVn)

loads <- cbind(loadsm, vips)
loads <- loads[order(loads$`model@vipVn`), ]
loads$features <- row.names(loads)
loads$features <- factor(loads$features, levels = loads$features)

# calculate enrichment----
pls.data <- data.frame(cbind(X_sel, lapply(y, as.character)))
colnames(pls.data)[dim(pls.data)[2]] <- "y"

pls.data.long <- tidyr::pivot_longer(pls.data, cols=colnames(X_sel), values_to = "value")
pls.data.long$value <- unlist(pls.data.long$value)
centroids <- pls.data.long %>% group_by(y, name) %>% summarise(val=mean(value))

enrich <- centroids %>% group_by(name) %>% summarise(group=y[which(val == max(val))])
enrich <- separate(enrich, col="name", into=c("assay", "antigen"), sep="_", remove = FALSE)
# enrich$antigen.num <- as.numeric(gsub("Pn", "", enrich$antigen))
enrich$antigen.num <- as.numeric(gsub("[^0-9.-]", "", enrich$antigen))
enrich$antigen.alph <- gsub("[0-9]+", "", enrich$antigen)

loads <- merge(loads, enrich, by.x = "features", by.y = "name")
loads <- loads[order(loads$antigen.num, loads$antigen.alph, loads$assay),]
loads <- as.data.frame(lapply(loads, unlist))

p.loads1 <- ggplot(data = loads, aes(x = features, y = p1)) +
  geom_bar(stat = "identity", aes(fill = group), colour = "black") +
  coord_flip() +
  scale_fill_manual(breaks = names(my_colors$group), values = as.vector(my_colors$group)) +
  scale_x_discrete(limits = rev(loads$features)) +
  theme_classic() +
  ylab("LV1 loadings") +
  xlab("") +
  theme(legend.position = "top")
print(p.loads1)

pdf(paste(saved.dir, paste(study, visit,  "pls-da_loading1.pdf", sep = "_"), sep = "/"), width = 4, height = 3)
print(p.loads1)
dev.off()

p.loads2 <- ggplot(data = loads, aes(x = features, y = p2)) +
  geom_bar(stat = "identity", aes(fill = group), colour = "black") +
  coord_flip() +
  scale_fill_manual(breaks = names(my_colors$group), values = my_colors$group) +
  scale_x_discrete(limits = rev(loads$features)) +
  theme_classic() +
  ylab("LV2 loadings") +
  xlab("") +
  theme(legend.position = "top")
print(p.loads2)

pdf(paste(saved.dir, paste(study, visit,  "pls-da_loading2.pdf", sep = "_"), sep = "/"), width = 4, height = 3)
print(p.loads2)
dev.off()

# networks-------
library(network)
library(stringr)
library(Hmisc)

temp <- as.data.frame(cbind(X, lapply(y, as.character)))
# assign name for Y
colnames(temp)[ncol(temp)] <- "vaccine"

q.cut <- 0.01
rho.cut <- .8

for (subtype in as.vector(unlist(unique(temp[ncol(temp)])))) {
  pred <- as.matrix(temp[temp$vaccine == subtype, 1:(ncol(temp) - 1)])

  res <- rcorr(as.matrix(pred), type = "spearman")

  resr <- as.data.frame(res$r)
  resr$selfeat <- row.names(resr)
  df.r <- subset(resr, resr$selfeat %in% sel_features)

  resp <- as.data.frame(res$P)
  resp$selfeat <- row.names(resp)
  df.p <- subset(resp, resr$selfeat %in% sel_features)

  dfr <- df.r %>% gather(feature, rho, -c(selfeat))
  dfp <- df.p %>% gather(feature, p, -c(selfeat))
  dfp$q <- p.adjust(dfp$p, method = "BH")

  identical(dfr$selfeat, dfp$selfeat)
  identical(dfr$feature, dfp$feature)

  df <- cbind(dfr, dfp[, c("p", "q")])
  # df$antigen <- as.vector(str_split_fixed(df$feature, "_", 2)[,2])

  df.sub <- subset(df, df$q < q.cut & abs(df$rho) > rho.cut)

  link.widths <- abs(df.sub$rho) * 5
  # link.widths <- ifelse(abs(df.sub$rho) > 0.5, 5, 1)

  link.clr <- ifelse(df.sub$rho > 0, "#e7d4e8", "#c7eae5")

  net <-
    network(
      df.sub,
      matrix.type = "edgelist",
      ignore.eval = FALSE,
      multiple = TRUE,
      directed = FALSE
    )

  node.cols <- network.vertex.names(net)
  fc.ind <- which(grepl("FcR|FcAR", node.cols))
  node.cols[fc.ind] <- "#1982c4"
  titer.ind <- which(grepl("Ig", node.cols))
  node.cols[titer.ind] <- "#ffffb3"

  rca.ind <- which(grepl("RCA", node.cols))
  node.cols[rca.ind] <- "#fb5607"
  sna.ind <- which(grepl("SNA", node.cols))
  node.cols[sna.ind] <- "#8ac926"

  func.ind <- which(grepl("ADCD|ADCP|ADNP|ADNK", node.cols))
  node.cols[func.ind] <- "#bebada"

  fname <-
    paste(saved.dir,
      paste(study, visit, subtype, "network.pdf", sep = "_"),
      sep = "/"
    )
  pdf(fname, width = 5, height = 5)
  print(
    plot.network(
      net,
      label = network.vertex.names(net),
      # vertex.col = "slategray",
      vertex.col = as.color(node.cols),
      # mode = "circle",
      edge.col = link.clr,
      edge.lwd = link.widths,
      pad = 2,
      label.cex = 0.5
    )
  )
  dev.off()
}

# cross validation------------------
y_pred <- predict_ropls(model, X_sel)
acc <- score_accuracy(y, y_pred)

print(paste("Performance on full data set:", round(acc, digits = 2), "accuracy"))

opts_sel <- list(threshold = 0.7, n_trials = 20, return_count = FALSE)
# update select function
opts_model <- list(n_LV = 2)

select_cv <- function(X, y) { return(select_repeat(X, y, selector = select_lasso, options = opts_sel)) }
train_cv <- function(X, y) { return(train_ropls(X, y, options = opts_model)) }

method <- list(select = select_cv,
              train = train_cv,
              predict = predict_ropls,
              score = score_accuracy)

opts <- list(n_folds = 5, rf_trials = 5, pt_trials = 5, n_trials = 5)
return.val <- validate(X, y, method, opts)
repeat.val <- validate_repeat(X, y, method, opts, n_trials = opts$n_trials)

# use return.val
accuracy <- c(return.val$cv_score, return.val$rf_scores, return.val$pt_scores)
model.method <- c(rep("model", 1), rep("random features", opts$n_trials), rep("permuted labels", opts$n_trials))
df.cv <- data.frame(accuracy, model.method)

# use repeat.val
accuracy <- c()
model.method <- c()

for (i in 1:opts$n_trials) {
  accuracy <- c(accuracy, repeat.val[[i]]$cv_score, repeat.val[[i]]$rf_scores, repeat.val[[i]]$pt_scores)
  model.method <- c(model.method, rep("model", 1), rep("random features", opts$rf_trials), rep("permuted labels", opts$pt_trials))
}
df.cv.repeat <- data.frame(accuracy, model.method)

my_comparisons = list(c("model", "random features"), c("model", "permuted labels"))
library(ggpubr)
fig.validation <- ggplot(data=df.cv.repeat, aes(x=model.method, y=accuracy)) +
  geom_violin( trim = FALSE) + ylim(0, 1) +
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1),
    geom = "pointrange", color = "black"
  ) +
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.95, 0.87)) +
  xlab("Model") +
  ylab("Accuracy") +
  theme_classic() + theme(
    legend.position = "top",
    axis.title = element_text(size = 10),
    axis.text.x = element_text(
      # angle = 45,
      color = "black",
      size = 10,
      vjust = 1.
    )
  ) 
print(fig.validation)
figurename <- paste(saved.dir, paste(study, visit, "pls_validation.pdf", sep = "_"), sep="/") 
ggsave(plot = fig.validation, dpi = 600, width = 4, height = 3,
       filename = figurename, device='pdf')
