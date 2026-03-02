source("../helper/common.R")

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

grid_options <- expand_grid(
  window = c(200, 250, 300, 350, 375, 400, 425, 450, 500, 600, 700, 800, 900, 1000),
  time_0 = c(F, T), time_end = c(F, T), during = T, frac = c(0, .5, .75)
) |>
  bind_rows(expand_grid(
    window = c(0),
    time_0 = c(F, T), time_end = c(F), during = T, frac = c(0)
  )) |>
  bind_rows(expand_grid(
    window = c(200, 250, 300, 350, 375, 400, 425, 450, 500, 600, 700, 800, 900, 1000),
    time_0 = c(T), time_end = c(T), during = T, frac = c(1)
  ))

rts <- rle_data %>%
  group_by(administration_id, trial_id, trial_order) %>%
  nest() |>
  cross_join(grid_options) |>
  mutate(rts = pmap(
    list(data, time_0, window, time_end, during, frac),
    \(d, t0, w, te, dur, fr){
      get_rt(d,
        t_0 = t0, window_length = w,
        t_end = te, window_mostly_region = dur, mostly_fraction = fr
      )
    }
  )) |>
  select(-data) |>
  unnest(cols = c(rts)) |>
  left_join(d_aoi %>%
    select(administration_id, dataset_name, subject_id, trial_id, trial_order, target_label) %>%
    distinct())

saveRDS(rts, "../cached_intermediates/3_rts.rds")
