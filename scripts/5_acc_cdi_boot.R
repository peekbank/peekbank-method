source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


downsample_acc_cdi <- function(t_start = -500, t_end = 4000, start_point, sample_down, iter) {
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
    cross_join(tibble(iteration = 1:iter)) |>
    group_by(administration_id, dataset_name, iteration) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name, iteration) |>
    summarize(mean_var = mean(accuracy, na.rm = T)) |>
    filter(!is.na(mean_var)) |>
    ungroup() |>
    left_join(cdi_data, by = c("administration_id", "dataset_name"))
}

params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)

accs_boot_cdi <- params |>
  mutate(iters = 1000) |>
  mutate(
    cdi = pmap(
      list(t_start, t_end, start_point, sample_down, iters),
      \(t_s, t_e, s_p, s_d, iters) {
        downsample_acc_cdi(t_s, t_e, s_p, s_d, iters) |>
          group_by(dataset_name, iteration) |>
          nest() |>
          mutate(corr = map2(data, dataset_name, \(d, dn) {
            calc_cdi(mutate(d, dataset_name = dn)) |> select(-any_of("dataset_name"))
          })) |>
          select(-data) |>
          unnest_wider(corr) |>
          ungroup() |>
          empirical_ci_cdi()
      }
    )
  ) |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "../cached_intermediates/5_acc_cdi_boot.rds")
