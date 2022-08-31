# library(showtext)
# 
# font_paths()  
# 
# trace(grDevices::png, exit = quote({
#   showtext::showtext_begin()
# }), print = FALSE)
# 
# 
# font_add("Arial", "../util/Arial.ttf")

theme_paperwhite <- function(
  base_size = 8L,
  base_family = "sans",
  face = c("plain", "bold"),
  aspect_ratio = NULL,
  legend_position = c("bottom", "right", "top", "none"),
  grid = FALSE,
  minimal = TRUE
) {
  
  face <- match.arg(face)
  legend_position <- match.arg(legend_position)
  
  
  gray <- "gray95"
  
  text <- element_text(
    family = "sans",
    face = face,
    colour = "#3d393d",
    size = 8L
  )
  
  # Include the grid lines.
  if (isTRUE(grid)) {
    panel_grid_major <- element_line(colour = gray, size = 0.5)
  } else {
    panel_grid_major <- element_blank()
  }
  
  # Remove panel border and axis ticks.
  if (isTRUE(minimal)) {
    axis_ticks <- element_line(colour = "#3d393d")
    panel_border <- element_blank()
    axis.line = element_line()
  } else {
    axis_ticks <- element_line(colour = "#3d393d")
    panel_border <- element_rect(colour = "#3d393d", fill = NA)
  }
  
  theme_linedraw(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      text = text,
      aspect.ratio = aspect_ratio,
      axis.line = element_line(),
      axis.text = text,
      # axis.text.x = element_text(angle = 90L, hjust = 1L, vjust = 0.5),
      axis.ticks = axis_ticks,
      panel.background = element_blank(),
      panel.border = panel_border,
      panel.grid.major = panel_grid_major,
      panel.grid.minor = element_blank(),
      legend.background = element_blank(),
      legend.position = legend_position,
      strip.background = element_rect(colour = NA, fill = "white"),
      strip.text = text,
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, face = "plain"),
      complete = TRUE,
      validate = TRUE
    )
}

theme_set(theme_paperwhite())
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = c("#424242","#a8404c","#024059","#71969F","#F2D6A2","#8B1D3B","#DD1822"))
}

scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = c("#424242","#a8404c","#024059","#71969F","#F2D6A2","#8B1D3B","#DD1822"))
}


col_pal <- c(colorRampPalette(colors = c("#424242", "#ffffff"))(49),
             "#ffffff","#ffffff",
             colorRampPalette(colors = c("#ffffff", "#a8404c"))(49))


col_pal1 <- c("#3d303d",
              "#443543",
              "#4a3a49",
              "#523d50",
              "#5b3b57",
              "#65345c",
              "#712860",
              "#7f195f",
              "#8e0b58",
              "#a00e50",
              "#b31241",
              "#c7162b")

col_pal2 <- c("#3d393d",
              "#443c43",
              "#4a3e49",
              "#523f50",
              "#5b3e57",
              "#653c5e",
              "#713864",
              "#7f3267",
              "#8e2c65",
              "#a0255d",
              "#b31e4a",
              "#c7162b")

col_pal3 <- c("#ffffff",
              "#fff1f8",
              "#ffe3ee",
              "#ffd0df",
              "#ffbacc",
              "#ffa1b5",
              "#ff869a",
              "#fc6b7e",
              "#f25264",
              "#e53b4d",
              "#d6273a",
              "#c7162b")

col_pal4 <- c("#ffffff",
              "#fff5f0",
              "#ffeae1",
              "#ffdbcd",
              "#ffc8b7",
              "#fdb29d",
              "#f99881",
              "#f37d66",
              "#eb634d",
              "#e14a37",
              "#d63424",
              "#ca2015")

col_pal5 <- c("#292222",
              "#2b2122",
              "#2d2021",
              "#322022",
              "#392125",
              "#432228",
              "#50232c",
              "#612430",
              "#772331",
              "#90212e",
              "#ad1d21",
              "#ca2015")

col_pal6 <- c("#707070",
              "#7a6d6a",
              "#826963",
              "#8a665d",
              "#926257",
              "#995e50",
              "#9f594a",
              "#a55544",
              "#ab503d",
              "#b14a37",
              "#b64431",
              "#bb3e2a",
              "#c03623",
              "#c52c1c",
              "#ca2015")

col_pal7 <- c("#fef2f1",
              "#fee4e2",
              "#fed6d2",
              "#fdc8c3",
              "#fcbab3",
              "#f9aca3",
              "#f69f94",
              "#f39184",
              "#ef8374",
              "#ea7465",
              "#e56655",
              "#df5746",
              "#d84736",
              "#d13626",
              "#ca2015")

col_pal8 <- c("#f7e4ea",
              "#f6ced7",
              "#f4b8c1",
              "#f19fa7",
              "#ed878b",
              "#e97071",
              "#e35a5a",
              "#dc4545",
              "#d33136",
              "#c81f2b",
              "#bd0f23")

col_pal9 <- c("#fcfbfb",
              "#fde6e8",
              "#fdd1d2",
              "#fbbebb",
              "#f9aca5",
              "#f69a8e",
              "#f08677",
              "#ea6e60",
              "#e1534a",
              "#d63434",
              "#ca212b",
              "#bd0f23")

col_pal10 <- colorRampPalette(colors = c("#43071E",
                                         "#691D32",
                                         "#923346",
                                         "#BD4B5C",
                                         "#D17486",
                                         "#E19EB0",
                                         "#F0C5D8",
                                         "#F8F0FE",
                                         "#C8D0EF",
                                         "#98B1DA",
                                         "#6A93C6",
                                         "#4272AE",
                                         "#31508C",
                                         "#1E356C",
                                         "#0E1949"))(100)

col_pal10 <- rev(col_pal10)


col_pal11 <- c("#f09648",
               "#ee9e5c",
               "#eba66f",
               "#e7ad81",
               "#e2b594",
               "#dcbda6",
               "#d5c4b9",
               "#cccccc",
               "#c2c5cf",
               "#b7bed2",
               "#abb7d5",
               "#a0b1d8",
               "#93aadb",
               "#86a3dd",
               "#779de0")
col_pal11 <- colorRampPalette(col_pal11)(100)
col_pal11 <- rev(col_pal11)