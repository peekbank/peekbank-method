library(tidyverse)
library(here)
library(tictoc)
library(ggthemes)
library(viridis)
library(metafor)
library(ggh4x)
knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed, cache = TRUE,
  message = FALSE, warning = FALSE, error = FALSE, echo = F,
  fig.height = 4, fig.width = 6
)
options(dplyr.summarise.inform = FALSE)


theme_set(theme_bw(base_size = 10))
theme_update(
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_line(),
  panel.grid.major.y = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  # panel.border = element_blank(),
  axis.line = element_line(),
  strip.text = element_text(face = "bold")
)

safe_rma <- function(d, col_est, col_upper, col_lower) {
  d <- d |>
    mutate(
      stdev = (.data[[col_upper]] - .data[[col_lower]]) / (1.96 * 2),
      var = stdev**2
    ) |>
    filter(!is.na(.data[[col_est]]), !is.na(var))
  if (nrow(d) > 2) {
    tryCatch(
      rma(d[[col_est]], d$var, control = list(maxiter = 1000)) |> summary() |> coef() |> select(estimate, ci.lb, ci.ub),
      error = function(e) tibble(estimate = NA, ci.lb = NA, ci.ub = NA)
    )
  } else {
    tibble(estimate = NA, ci.lb = NA, ci.ub = NA)
  }
}

do_cdi_meta <- function(dataset, grouping_factors) {
  dataset |>
    group_by(across(all_of(grouping_factors))) |>
    nest() |>
    mutate(
      comp = map(data, \(d){
        safe_rma(d, "comp_est", "comp_upper", "comp_lower")
      }),
      prod = map(data, \(d){
        safe_rma(d, "prod_est", "prod_upper", "prod_lower")
      }),
      age = map(data, \(d){
        safe_rma(d, "age_est", "age_upper", "age_lower")
      }),
    ) |>
    select(-data) |>
    unnest(c(comp, prod, age), names_sep = "_")
}

do_meta <- function(dataset, grouping_factors) {
  dataset |>
    group_by(across(all_of(grouping_factors))) |>
    nest() |>
    mutate(corr = map(data, \(d){
      safe_rma(d, "est", "upper", "lower")
    })) |>
    unnest(corr) |>
    select(-data)
}

regroup_rt <- function(df) {
  df |>
    filter(time_0, time_end, during, frac == 1) |>
    mutate(
      type = case_when(
        str_detect(measure, "first_launch_rt") ~ "first_launch",
        str_detect(measure, "last_launch_rt") ~ "last_launch",
        str_detect(measure, "land_rt") ~ "land"
      ),
      logged = case_when(
        str_detect(measure, "log") ~ "log",
        T ~ "raw"
      ),
      trimming = case_when(
        str_detect(measure, "trim_first") ~ "trim_first",
        str_detect(measure, "trim_last") ~ "trim_last",
        T ~ "untrimmed"
      )
    ) |>
    filter(trimming != "trim_last")
}

name_rt <- function(df) {
  df |>
    mutate(approach = case_when(
      type == "land" & trimming == "untrimmed" & logged == "log" ~ "Log Land RT",
      type == "land" & trimming == "untrimmed" & logged == "raw" ~ "Raw Land RT",
      type == "first_launch" & trimming == "trim_first" & logged == "log" ~ "Log Launch RT",
      type == "first_launch" & trimming == "trim_first" & logged == "raw" ~ "Raw Launch RT",
    )) |>
    filter(!is.na(approach))
}

combined_metrics <- function(icc, cdi, trt) {
  icc |>
    mutate(Type = "ICC reliability") |>
    bind_rows(cdi |> rename(estimate = prod_estimate, ci.lb = prod_ci.lb, ci.ub = prod_ci.ub) |> mutate(Type = "Corr. with CDI Prod.")) |>
    bind_rows(cdi |> rename(estimate = comp_estimate, ci.lb = comp_ci.lb, ci.ub = comp_ci.ub) |> mutate(Type = "Corr. with CDI Comp.")) |>
    bind_rows(trt |> mutate(Type = "Test-retest reliability"))
}

combined_metrics_age <- function(icc, cdi) {
  icc |>
    mutate(Type = "ICC reliability") |>
    bind_rows(cdi |> rename(estimate = prod_estimate, ci.lb = prod_ci.lb, ci.ub = prod_ci.ub) |> mutate(Type = "Corr. with CDI Prod.")) |>
    bind_rows(cdi |> rename(estimate = comp_estimate, ci.lb = comp_ci.lb, ci.ub = comp_ci.ub) |> mutate(Type = "Corr. with CDI Comp."))
}

order_age <- function(df) {
  df |> mutate(age_bin = factor(age_bin, levels = c("<18", "18-24", "24-36", ">=36")))
}
