#' Plot correlation as a chord diagram for given columns A and columns B
#'
#' @param df: n_samples x n_features matrix
#' @param feat.cols: features group A, vector of strings
#' @param title: features group B, vector of strings
#'
#' @return none
#' @export
#'
draw_cor_chord <- function(df, feat.cols, title, pcut = 0.1, rcut = 0.5, colors = NULL) {
  if (is.null(colors)) {
    colors <- list(
      max.color = "#40004b",
      half.max.color = "#e7d4e8",
      half.min.color = "#c7eae5",
      min.color = "#003c30"
    )
  }
  # browser()
  circ <-
    subset(df, abs(df$cor) > rcut & df$p < pcut)

  circ <- circ[, c("from", "to", "cor")]

  colnames(circ) <- c("from", "to", "value")
  circ$from <- as.factor(circ$from)
  circ$to <- as.factor(circ$to)
  circ$value <- as.numeric(circ$value)
  sections <- c(levels(circ$to), levels(circ$from))

  col_fun <- circlize::colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c(
      colors$min.color,
      colors$half.min.color,
      "white",
      colors$half.max.color,
      colors$max.color
    )
  )

  pdf(title, width = 3, height = 2.5)
  circlize::chordDiagram(
    circ,
    order = sections,
    grid.col = feat.cols,
    col = col_fun,
    annotationTrack = c("grid"),
    preAllocateTracks = list(track.height = max(strwidth(
      unlist(dimnames(circ))
    )))
  )

  circlize::circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim <- circlize::get.cell.meta.data("xlim")
      xplot <- circlize::get.cell.meta.data("xplot")
      ylim <- circlize::get.cell.meta.data("ylim")
      sector.name <- circlize::get.cell.meta.data("sector.index")
      if (abs(xplot[2] - xplot[1]) < 5.5) {
        circlize::circos.text(
          mean(xlim),
          ylim[1],
          sector.name,
          facing = "clockwise",
          cex = 0.18,
          niceFacing = TRUE,
          adj = c(0, 0.5),
          col = "black"
        )
      } else {
        circlize::circos.text(
          mean(xlim),
          ylim[1],
          sector.name,
          facing = "clockwise",
          cex = 0.22,
          niceFacing = TRUE,
          adj = c(0, 0.5),
          col = "black"
        )
      }
    },
    bg.border = NA
  )
  title(main = list(
    paste0(""),
    cex = 0.8,
    font = 1,
    col = "black"
  ))
  dev.off()
}

#' Generate palette for a chord diagram
#'
#' @param col.feat: features with formatted names
#'
#' @return none
#' @export
#'
gen_sec_color <- function(col.feat, color_col = "assay", dlim = "_") {
  df <- as.data.frame(col.feat)
  df <-
    separate(
      df,
      "col.feat",
      into = c("assay", "antigen"),
      sep = dlim,
      remove = FALSE
    )

  my.pal <-
    data.frame(
      color = circlize::rand_color(length(unique(df[[color_col]])), transparency = 0.1),
      assay = unique(df[[color_col]])
    )
  temp <- merge(df, my.pal, by = color_col, all.x = TRUE)
  feat.cols <- structure(temp$color, names = temp$col.feat)
  return(feat.cols)
}
