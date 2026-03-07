source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


downsample_acc_cdi <- function(t_start = -500, t_end = 4000, start_point, sample_down) {
  d_aoi |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name) |>
    mutate(count = n()) |>
    filter(count >= start_point) |>
    select(-count) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name) |>
    summarize(mean_var = mean(accuracy, na.rm = T)) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data)
}

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "do_cdi", "boot_cdi")
)



params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)

accs_boot_cdi <- params |>
  mutate(summary_data = pmap(list(t_start, t_end, start_point, sample_down), \(t_s, t_e, s_p, s_d) downsample_acc_cdi(t_s, t_e, s_p, s_d))) |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, boot_cdi)) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "../cached_intermediates/5_acc_cdi_boot.rds")
