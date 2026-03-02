source("cl_helper.R")

d_aoi <- readRDS("d_aoi.Rds")

cdi_data <- readRDS("cdi.rds")



acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 4000),
)


do_cdi_accuracy <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    group_by(administration_id) |>
    summarize(mean_var = mean(accuracy, na.rm = T)) |>
    left_join(cdi_data) |>
    ungroup() |>
    summarise(
      cor_comp = ifelse(sum(!is.na(comp) & !is.na(mean_var)) > 2, cor.test(mean_var, comp)$estimate, as.numeric(NA)),
      cor_prod = ifelse(sum(!is.na(prod) & !is.na(mean_var)) > 2, cor.test(mean_var, prod)$estimate, as.numeric(NA)),
      cor_age = ifelse(sum(!is.na(age) & !is.na(mean_var)) > 2, cor.test(mean_var, age)$estimate, as.numeric(NA))
    )
  cors <- c(summ$cor_comp[1], summ$cor_prod[1], summ$cor_age[1])
  names(cors) <- c("cor_comp", "cor_prod", "cor_age")

  return(cors)
}

acc_cdi <- function(t_start = -500, t_end = 4000) {
  df <- d_aoi |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    group_by(dataset_name, dataset_id, administration_id, target_label) |>
    filter(!is.na(accuracy))

  # compute ICCs
  df |>
    group_by(dataset_name) |>
    nest() |>
    # head(5) |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi_accuracy, 2000)
      if (is.na(b$t0[1])) {
        comp_lower <- NA
        comp_upper <- NA
      } else {
        ci_comp <- boot::boot.ci(b, index = 1, type = "basic")
        comp_lower <- ci_comp$basic[4]
        comp_upper <- ci_comp$basic[5]
      }
      if (is.na(b$t0[2])) {
        prod_lower <- NA
        prod_upper <- NA
      } else {
        ci_prod <- boot::boot.ci(b, index = 2, type = "basic")
        prod_lower <- ci_prod$basic[4]
        prod_upper <- ci_prod$basic[5]
      }
      if (is.na(b$t0[3])) {
        age_lower <- NA
        age_upper <- NA
      } else {
        ci_age <- boot::boot.ci(b, index = 3, type = "basic")
        age_lower <- ci_age$basic[4]
        age_upper <- ci_age$basic[5]
      }
      tibble(
        comp_est = b$t0[1], comp_lower = comp_lower, comp_upper = comp_upper,
        prod_est = b$t0[2], prod_lower = prod_lower, prod_upper = prod_upper,
        age_est = b$t0[3], age_lower = age_lower, age_upper = age_upper,
      )
    })) |>
    select(-data) |>
    unnest(corr)
}


cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_copy(cluster, "do_cdi_accuracy")
cluster_copy(cluster, "d_aoi")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "acc_cdi")

accs_boot_cdi <- acc_params |>
  partition(cluster) |>
  # head(1) |>
  mutate(cdi = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi(t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "3_accs_boot_cdi.rds")


do_cdi_bc_acc <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    group_by(administration_id) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T)) |>
    left_join(cdi_data) |>
    ungroup() |>
    summarise(
      cor_comp = ifelse(sum(!is.na(comp) & !is.na(mean_var)) > 2, cor.test(mean_var, comp)$estimate, as.numeric(NA)),
      cor_prod = ifelse(sum(!is.na(prod) & !is.na(mean_var)) > 2, cor.test(mean_var, prod)$estimate, as.numeric(NA)),
      cor_age = ifelse(sum(!is.na(age) & !is.na(mean_var)) > 2, cor.test(mean_var, age)$estimate, as.numeric(NA))
    )
  cors <- c(summ$cor_comp[1], summ$cor_prod[1], summ$cor_age[1])
  names(cors) <- c("cor_comp", "cor_prod", "cor_age")

  return(cors)
}

bc_acc_cdi <- function(b_start, b_end, t_start, t_end) {
  # get trials with some baseline
  baseline_lengths <- d_aoi |>
    group_by(dataset_name, trial_id) |>
    summarise(t_min = min(t_norm))

  # get baseline corrected accuracies for all trials with ANY baseline info
  df <- d_aoi |>
    left_join(baseline_lengths) |>
    filter(t_min < 0) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    filter(!is.na(bc_accuracy)) |>
    group_by(dataset_id, dataset_name, administration_id, target_label)

  # compute ICCs
  df |>
    group_by(dataset_name) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi_bc_acc, 2000)
      if (is.na(b$t0[1])) {
        comp_lower <- NA
        comp_upper <- NA
      } else {
        ci_comp <- boot::boot.ci(b, index = 1, type = "basic")
        comp_lower <- ci_comp$basic[4]
        comp_upper <- ci_comp$basic[5]
      }
      if (is.na(b$t0[2])) {
        prod_lower <- NA
        prod_upper <- NA
      } else {
        ci_prod <- boot::boot.ci(b, index = 2, type = "basic")
        prod_lower <- ci_prod$basic[4]
        prod_upper <- ci_prod$basic[5]
      }
      if (is.na(b$t0[3])) {
        age_lower <- NA
        age_upper <- NA
      } else {
        ci_age <- boot::boot.ci(b, index = 3, type = "basic")
        age_lower <- ci_age$basic[4]
        age_upper <- ci_age$basic[5]
      }
      tibble(
        comp_est = b$t0[1], comp_lower = comp_lower, comp_upper = comp_upper,
        prod_est = b$t0[2], prod_lower = prod_lower, prod_upper = prod_upper,
        age_est = b$t0[3], age_lower = age_lower, age_upper = age_upper,
      )
    })) |>
    select(-data) |>
    unnest(corr)
}


cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_copy(cluster, "do_cdi_bc_acc")
cluster_copy(cluster, "d_aoi")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "bc_acc_cdi")

bc_acc_params <- expand_grid(
  t_start = 500,
  t_end = c(4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

bc_accs_boot_cdi <- bc_acc_params |>
  partition(cluster) |>
  mutate(cdi = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_acc_cdi(b_s, b_e, t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(bc_accs_boot_cdi, "7_bc_accs_boot_cdi.rds")


do_cdi_rt <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    group_by(administration_id) |>
    summarize(rt = mean(rt, na.rm = T)) |>
    left_join(cdi_data) |>
    ungroup() |>
    summarise(
      cor_comp = ifelse(sum(!is.na(comp) & !is.na(rt)) > 2, cor.test(rt, comp)$estimate, as.numeric(NA)),
      cor_prod = ifelse(sum(!is.na(prod) & !is.na(rt)) > 2, cor.test(rt, prod)$estimate, as.numeric(NA)),
      cor_age = ifelse(sum(!is.na(age) & !is.na(rt)) > 2, cor.test(rt, age)$estimate, as.numeric(NA))
    )
  cors <- c(summ$cor_comp[1], summ$cor_prod[1], summ$cor_age[1])
  names(cors) <- c("cor_comp", "cor_prod", "cor_age")

  return(cors)
}



cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
# cluster_library(cluster, "agreement")
# cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "do_cdi_rt")
cluster_copy(cluster, "cdi_data")

rts <- readRDS("rts.rds")


d_rt_dt <- rts |>
  filter(window %in% c(200, 350, 375, 400, 425, 450, 500)) |>
  filter(frac %in% c(0, 1)) |>
  filter(time_0 == T, time_end == T) |>
  filter(shift_type == "D-T") |>
  mutate(
    land_rt = rt,
    first_launch_rt = shift_start_rt,
    # last_launch_rt = last_shift_rt
  ) |>
  mutate(
    across(c("land_rt", "first_launch_rt"), log, .names = "log_{.col}"),
    across(
      c(
        "first_launch_rt",
        "log_first_launch_rt",
      ),
      ~ ifelse(shift_length <= 600, .x, NA),
      .names = "trim_first_{.col}"
    ),
  ) |>
  select(-rt, -shift_start_rt, -last_shift_rt) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label)

d_rt_dt_long <- d_rt_dt |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt") |>
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
      T ~ "untrimmed"
    )
  ) |>
  mutate(approach = case_when(
    type == "land" & trimming == "untrimmed" & frac == "0" ~ "generous",
    type == "land" & trimming == "untrimmed" & frac == "1" ~ "trad_land",
    type == "first_launch" & trimming == "trim_first" & frac == "1" ~ "trad_launch"
  )) |>
  filter(!is.na(approach))

bootstrap_cdi <- d_rt_dt_long |>
  ungroup() |>
  select(administration_id, rt, approach, window, logged, dataset_name) |>
  group_by(approach, window, logged, dataset_name) |>
  nest() |>
  partition(cluster) |>
  mutate(corr = map(data, \(d) {
    b <- boot::boot(d, do_cdi_rt, 2000)
    if (is.na(b$t0[1])) {
      comp_lower <- NA
      comp_upper <- NA
    } else {
      ci_comp <- boot::boot.ci(b, index = 1, type = "basic")
      comp_lower <- ci_comp$basic[4]
      comp_upper <- ci_comp$basic[5]
    }
    if (is.na(b$t0[2])) {
      prod_lower <- NA
      prod_upper <- NA
    } else {
      ci_prod <- boot::boot.ci(b, index = 2, type = "basic")
      prod_lower <- ci_prod$basic[4]
      prod_upper <- ci_prod$basic[5]
    }
    if (is.na(b$t0[3])) {
      age_lower <- NA
      age_upper <- NA
    } else {
      ci_age <- boot::boot.ci(b, index = 3, type = "basic")
      age_lower <- ci_age$basic[4]
      age_upper <- ci_age$basic[5]
    }
    tibble(
      comp_est = b$t0[1], comp_lower = comp_lower, comp_upper = comp_upper,
      prod_est = b$t0[2], prod_lower = prod_lower, prod_upper = prod_upper,
      age_est = b$t0[3], age_lower = age_lower, age_upper = age_upper,
    )
  })) |>
  collect() |>
  select(-data) |>
  unnest(corr)

saveRDS(bootstrap_cdi, "4_rt_cdi_boot.rds")
