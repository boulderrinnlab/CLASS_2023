theme_paperwhite <- function(
  base_size = 14L,
  base_family = "",
  face = c("bold", "plain"),
  aspect_ratio = NULL,
  legend_position = c("right", "bottom", "top", "none"),
  grid = FALSE,
  minimal = FALSE
) {
  
  face <- match.arg(face)
  legend_position <- match.arg(legend_position)
  
  
  gray <- "gray95"
  
  text <- element_text(
    family = base_family,
    face = face,
    colour = "#424242"
  )
  
  # Include the grid lines.
  if (isTRUE(grid)) {
    panel_grid_major <- element_line(colour = gray, size = 0.5)
  } else {
    panel_grid_major <- element_blank()
  }
  
  # Remove panel border and axis ticks.
  if (isTRUE(minimal)) {
    axis_ticks <- element_blank()
    panel_border <- element_blank()
  } else {
    axis_ticks <- element_line(colour = "#424242")
    panel_border <- element_rect(colour = "#424242", fill = NA)
  }
  
  theme_linedraw(
    base_size = base_size,
    base_family = base_family
  ) +
    theme(
      text = text,
      aspect.ratio = aspect_ratio,
      axis.line = element_blank(),
      axis.text = text,
      axis.text.x = element_text(angle = 90L, hjust = 1L, vjust = 0.5),
      axis.ticks = axis_ticks,
      panel.background = element_blank(),
      panel.border = panel_border,
      panel.grid.major = panel_grid_major,
      panel.grid.minor = element_blank(),
      legend.background = element_blank(),
      legend.position = legend_position,
      strip.background = element_rect(colour = NA, fill = "white"),
      strip.text = text,
      complete = TRUE,
      validate = TRUE
    )
}

theme_set(
  theme_paperwhite(
    base_size = 14L,
    legend_position = "right"
  )
)