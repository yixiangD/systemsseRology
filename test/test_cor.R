# read data
data <- readxl::read_excel("data/data_truncated.xlsx")
# compute correlation
col.sel <- colnames(data)[grepl("AUC", colnames(data))]

library(dplyr)
source("R/compute_correlation.R")
df <- get_cor_tri(data, col.sel)
df <- get_cor_col(data, col.sel, col.sel)
df$signif <- systemsseRology::pval_to_asterisk(df$p)
print(dim(df))
library(ggplot2)
fig.heatmap <-
    ggplot(data = df, aes(x = from, y = to, fill = cor)) +
    geom_tile(
      width = 1,
      height = 1,
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(aes(label = signif), size = 3) +
    scale_y_discrete(limits = unique(df$from)) +
    scale_x_discrete(limits = unique(df$to)) +
    labs(x = "From", y = "To") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 10),
      axis.text.y = element_text(
        size = 10,
        color = "black",
        hjust = 1,
        vjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 45,
        color = "black",
        size = 10,
        hjust = 1,
        vjust = 1.0
      ),
      legend.position = "top"
    ) +
    theme(strip.background = element_blank()) +
    scale_fill_gradient2(
      low = "#6c757d",
      high = "#f94144",
      mid = "white",
      na.value = "black",
      limits = c(-1., 1.),
      midpoint = 0,
      name = paste("", sep = "")
    )
print(fig.heatmap)
