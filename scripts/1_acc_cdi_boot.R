source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


acc_cdi <- function(t_start = -500, t_end = 4000) {
  d_aoi |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    group_by(dataset_name, dataset_id, administration_id, target_label) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name) |>
    summarize(mean_var = mean(accuracy, na.rm = T)) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data, by = c("dataset_name", "administration_id"))
}



acc_cdi_age <- function(t_start = -500, t_end = 4000) {
  d_aoi_age |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, age_bin, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, age_bin) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, age_bin) |>
    summarize(mean_var = mean(accuracy, na.rm = T)) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data, by = c("dataset_name", "administration_id"))
}


acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)

accs_cdi <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi(t_s, t_e))) |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi, "../cached_intermediates/1_acc_cdi.rds")

accs_cdi_byage <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi_age(t_s, t_e))) |>
  mutate(cdi = map(summary_data, \(d) calc_cdi(d, by_age = TRUE))) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi_byage, "../cached_intermediates/1_acc_cdi_byage.rds")
