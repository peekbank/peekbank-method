source("../helper/common.R")
source("../helper/rt_helper.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> left_join(rt_params)

# rts_weird <- readRDS("../cached_intermediates/3_rts.rds") |> filter(min_rt == 400)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

d_rt_dt <- preprocess_rt_dt(rts) |>
  group_by(dataset_name, time_0, min_rt, max_rt, time_end, during, frac, administration_id, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var)) |>
  left_join(cdi_data, by = c("dataset_name", "administration_id"))

d_rt_dt_byage <- preprocess_rt_dt(rts) |>
  inner_join(age_bin_cutoff) |>
  group_by(dataset_name, time_0, min_rt, max_rt, time_end, during, frac, administration_id, age_bin, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var)) |>
  left_join(cdi_data, by = c("dataset_name", "administration_id"))

# d_rt_dt_weird <- preprocess_rt_dt(rts_weird) |>
#   filter(measure %in% c("log_land_rt", "land_rt")) |>
#   group_by(dataset_name, time_0, min_rt, max_rt, time_end, during, frac, administration_id, measure) |>
#   summarize(mean_var = mean(rt, na.rm = T)) |>
#   filter(!is.na(mean_var)) |>
#   left_join(cdi_data)

rt_cdi <- d_rt_dt |>
  group_by(time_0, min_rt, max_rt, time_end, during, frac, measure) |>
  nest() |>
  mutate(cdi = map(data, calc_cdi)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_cdi, "../cached_intermediates/3_rt_cdi.rds")

rt_cdi_byage <- d_rt_dt_byage |>
  group_by(age_bin, time_0, min_rt, max_rt, time_end, during, frac, measure) |>
  nest() |>
  mutate(cdi = map(data, calc_cdi)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_cdi_byage, "../cached_intermediates/3_rt_cdi_byage.rds")

# rt_cdi_weird <- d_rt_dt_weird |>
#   group_by(time_0, min_rt, max_rt, time_end, during, frac, measure) |>
#   nest() |>
#   mutate(cdi = map(data, calc_cdi)) |>
#   select(-data) |>
#   unnest(cdi)
#
# saveRDS(rt_cdi_weird, "../cached_intermediates/3_rt_cdi_weird.rds")
