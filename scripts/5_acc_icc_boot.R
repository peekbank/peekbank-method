source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

# Summarize timepoint-level data to trial-level accuracy.
# Includes age_bin in grouping if present in the data.

summarized_accuracy <- d_aoi |>
  filter(t_norm > 400, t_norm < 3000) |>
  group_by(across(all_of(c(
    "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
  ))), across(any_of("age_bin"))) |>
  summarise(
    accuracy = mean(correct, na.rm = TRUE),
    prop_data = mean(!is.na(correct)),
    .groups = "drop"
  ) |>
  filter(!is.na(accuracy)) |>
  group_by(across(all_of(c(
    "dataset_name", "dataset_id", "administration_id"
  ))), across(any_of("age_bin"))) |>
  mutate(count = n())

rm(d_aoi)

downsample_summarize_accuracy <- function(start_point, sample_down, iter) {
  summarized_accuracy |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iter)) |>
    group_by(iteration, across(all_of(c(
      "dataset_name", "dataset_id", "administration_id"
    ))), across(any_of("age_bin"))) |>
    slice_sample(n = sample_down) |>
    group_by(iteration, across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup() |>
    summarize_icc_resamples("accuracy")
}

params <- expand_grid(
  start_point = c(5, 10, 15, 20),
  sample_down = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)


accs_boot <- params |>
  mutate(iters = 500) |>
  mutate(corr = pmap(list(start_point, sample_down, iters), \(s_p, s_d, iters) downsample_summarize_accuracy(s_p, s_d, iters))) |>
  unnest(corr)

saveRDS(accs_boot, "../cached_intermediates/5_acc_icc_boot.rds")
