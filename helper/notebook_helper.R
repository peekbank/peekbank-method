library(tidyverse)
library(here)
library(tictoc)
library(ggthemes)
library(viridis)
library(metafor)
knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed, cache = TRUE,
  message = FALSE, warning = FALSE, error = FALSE
)
options(dplyr.summarise.inform = FALSE)


theme_set(theme_bw(base_size = 10))
theme_update(
  panel.grid = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(),
  strip.text = element_text(face = "bold")
)

do_cdi_meta <- function(dataset, grouping_factors) {
  dataset |>
    mutate(
      comp_stdev = (comp_upper - comp_lower) / (1.96 * 2),
      comp_var = comp_stdev**2,
      prod_stdev = (prod_upper - prod_lower) / (1.96 * 2),
      prod_var = prod_stdev**2,
      age_stdev = (age_upper - age_lower) / (1.96 * 2),
      age_var = age_stdev**2
    ) |>
    group_by(across(all_of(grouping_factors))) |>
    nest() |>
    mutate(
      comp = map(data, \(d){
        rma(d$comp_est, d$comp_var) |>
          summary() |>
          coef()
      }),
      prod = map(data, \(d){
        rma(d$prod_est, d$prod_var) |>
          summary() |>
          coef()
      }),
      age = map(data, \(d){
        rma(d$age_est, d$age_var) |>
          summary() |>
          coef()
      })
    ) |>
    select(-data) |>
    unnest(c(comp, prod, age), names_sep = "_")
}

do_meta <- function(dataset, grouping_factors) {
  dataset |>
    filter(!is.na(lower), !is.na(upper)) |>
    mutate(
      stdev = (upper - lower) / (1.96 * 2),
      var = stdev**2,
    ) |>
    group_by(across(all_of(grouping_factors))) |>
    nest() |>
    mutate(corr = map(data, \(d){
      rma(d$est, d$var) |>
        summary() |>
        coef()
    })) |>
    unnest(corr)
}
