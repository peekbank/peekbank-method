knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed, cache = TRUE,
  message = FALSE, warning = FALSE, error = FALSE
)
options(dplyr.summarise.inform = FALSE)


.font <- "Source Sans Pro"
theme_set(theme_bw(base_size = 14, base_family = .font))
theme_set(theme_bw(base_size = 10))
theme_update(
  panel.grid = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(),
  strip.text = element_text(face = "bold")
)
