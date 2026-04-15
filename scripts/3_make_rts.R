source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

source("../helper/rt_helper.R")

rle_data <- d_aoi %>%
  filter(
    any(t_norm == 0), # must have data at 0
    t_norm >= 0
  ) %>% # only pass data after 0
  group_by(administration_id, trial_id, trial_order, dataset_name) %>%
  reframe(
    lengths = rle(aoi)$lengths,
    values = rle(aoi)$values
  )

grid_options <- rt_params |>
  select(-max_rt) |>
  distinct()

rts <- rle_data %>%
  group_by(administration_id, trial_id, trial_order) %>%
  nest() |>
  cross_join(grid_options) |>
  mutate(rts = pmap(
    list(data, time_0, min_rt, time_end, during, frac),
    \(d, t0, minr, te, dur, fr){
      get_rt(d,
        t_0 = t0, window_length_ms = minr,
        t_end = te, window_mostly_region = dur, mostly_fraction = fr
      )
    }
  )) |>
  select(-data) |>
  unnest(cols = c(rts)) |>
  left_join(d_aoi %>%
    select(administration_id, dataset_name, subject_id, trial_id, trial_order, target_label, age) %>%
    distinct())

saveRDS(rts, "../cached_intermediates/3_rts.rds")


weird_rt_grid <- rt_params_weird |>
  select(-max_rt) |>
  distinct()


weird_rts <-
  rle_data %>%
  group_by(administration_id, trial_id, trial_order) %>%
  nest() |>
  cross_join(weird_rt_grid) |>
  mutate(rts = pmap(
    list(data, time_0, min_rt, time_end, during, frac),
    \(d, t0, minr, te, dur, fr){
      get_rt(d,
        t_0 = t0, window_length_ms = minr,
        t_end = te, window_mostly_region = dur, mostly_fraction = fr
      )
    }
  )) |>
  select(-data) |>
  unnest(cols = c(rts)) |>
  left_join(d_aoi %>%
    select(administration_id, dataset_name, subject_id, trial_id, trial_order, target_label, age) %>%
    distinct())

saveRDS(weird_rts, "../cached_intermediates/3_weird_rts.rds")
