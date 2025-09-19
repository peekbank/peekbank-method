# source me for all analyses

library(here)
library(tidyverse)
library(multidplyr)
# library(ggpmisc)
library(ggrepel)
# library(ggthemes)
library(viridis)
# library(cowplot)
# remotes::install_github("jmgirard/agreement")
library(agreement)
library(tictoc)

# Seed for random number generation
set.seed(42)
knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed, cache = TRUE,
  message = FALSE, warning = FALSE, error = FALSE
)
options(dplyr.summarise.inform = FALSE)

# load data
# d_aoi <- readRDS(file = here("cached_intermediates","1_d_aoi.Rds"))

# key ICC function
get_icc <- function(x, column = "accuracy", object = "stimulus") {
  if (object == "stimulus") {
    iccs <- dim_icc(x,
      model = "2A",
      type = "consistency",
      unit = "average",
      object = target_label,
      rater = administration_id,
      trial = repetition,
      score = {{ column }},
      bootstrap = 0
    )
  } else {
    iccs <- dim_icc(x,
      model = "2A",
      type = "consistency",
      unit = "average",
      object = administration_id,
      rater = target_label,
      trial = repetition,
      score = {{ column }},
      bootstrap = 0
    )
  }

  return(iccs$Inter_ICC)
}
bootstrap_icc <- function(x, column = "accuracy", bootstrap = 2000) {
  iccs <- dim_icc(x,
    model = "2A",
    type = "consistency",
    unit = "average",
    object = administration_id,
    rater = target_label,
    trial = repetition,
    score = {{ column }},
    bootstrap = bootstrap
  )
  # View(iccs$boot_results[["t"]])
  t <- iccs$boot_results[["t"]][, 6]
  t_no_na <- t[!is.na(t)]
  # print(t_no_na)
  # View(iccs$boot_results[["t0"]])
  t0 <- iccs$boot_results[["t0"]][6][1]
  # print(t0)
  if (is.na(t0)) {
    return(tibble(est = NA, lower = NA, upper = NA))
  }
  if (length(t_no_na) == 0) {
    return(tibble(est = iccs$Inter_ICC, lower = NA, upper = NA))
  }
  iccs$boot_results$R <- length(t_no_na)
  ci <- boot::boot.ci(boot.out = iccs$boot_results, t = t_no_na, t0 = t0, type = "basic")
  return(tibble(est = iccs$Inter_ICC, lower = ci$basic[4], upper = ci$basic[5]))
}
options(dplyr.summarise.inform = FALSE)


# .font <- "Source Sans Pro"
# theme_set(theme_bw(base_size = 14, base_family = .font))
theme_set(theme_bw(base_size = 10))
theme_update(
  panel.grid = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(),
  strip.text = element_text(face = "bold")
)

# helper function for scaling to the variance of a particular age group
age_scale <- function(x, age, min = 16, max = 20) {
  (x - mean(x, na.rm = TRUE)) / sd(x[age >= min & age <= max], na.rm = TRUE)
}
