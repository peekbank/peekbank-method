source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

d_aoi_bc <- make_baseline_corrected(d_aoi)
d_aoi_bc_age <- make_baseline_corrected(d_aoi_age)

bc_acc_cdi <- function(b_start = -2000, b_end = 0,
                       t_start = 500, t_end = 4000) {
  d_aoi_bc |>
    group_by(dataset_name, dataset_id, administration_id, trial_id) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    group_by(administration_id, dataset_name) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T)) |>
    filter(!is.na(mean_var))
}



bc_acc_cdi_age <- function(b_start = -2000, b_end = 0,
                           t_start = 500, t_end = 4000) {
  d_aoi_bc_age |>
    group_by(dataset_name, dataset_id, administration_id, trial_id, age_bin) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    group_by(administration_id, dataset_name, age_bin) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T)) |>
    filter(!is.na(mean_var))
}

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "do_cdi", "d_aoi_bc", "d_aoi_bc_age", "cdi_data", "boot_cdi")
)


bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

bc_boot_cdi <- bc_acc_params |>
  mutate(summary_data = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_acc_cdi(b_s, b_e, t_s, t_e))) |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, boot_cdi)) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(bc_boot_cdi, "../cached_intermediates/2_bc_cdi_boot.rds")

bc_boot_cdi_age <- bc_acc_params |>
  mutate(summary_data = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_acc_cdi_age(b_s, b_e, t_s, t_e))) |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, \(d) boot_cdi(d, by_age = TRUE))) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(bc_boot_cdi_age, "../cached_intermediates/2_bc_cdi_boot_byage.rds")
