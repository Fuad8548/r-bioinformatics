# R/plotting_theme.R - Custom publication theme for ggplot2
library(ggplot2)

theme_bio <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = 'bold'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'bottom'
    )
}
