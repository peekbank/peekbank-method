# source me for all analyses

library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tibble)
library(agreement)
library(multidplyr)
# because we don't have full on tidyverse on the cluster, we specify the parts we actually need

# Seed for random number generation
set.seed(42)


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
