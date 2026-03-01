source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- d_aoi |>
  filter(!is.na(correct)) |>
  distinct(administration_id, age, dataset_name) |>
  mutate(age_bin = case_when(
    age < 18 ~ "<18",
    age < 24 ~ "18-24",
    age < 36 ~ "24-36",
    age >= 36 ~ ">=36"
  )) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= 5) |>
  ungroup()

d_aoi_age <- d_aoi |> inner_join(age_bin_cutoff)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

do_cdi <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
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

boot_cdi <- function(data) {
  data |>
    group_by(dataset_name) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi, 2000)
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
    boot_cdi()
}


boot_cdi_age <- function(data) {
  data |>
    group_by(dataset_name, age_bin) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi, 2000)
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
    boot_cdi_age()
}



cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_copy(cluster, "do_cdi")
cluster_copy(cluster, "d_aoi")
cluster_copy(cluster, "d_aoi_age")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "boot_cdi")
cluster_copy(cluster, "acc_cdi")
cluster_copy(cluster, "boot_cdi_age")
cluster_copy(cluster, "acc_cdi_age")


acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)

accs_boot_cdi <- acc_params |>
  partition(cluster) |>
  mutate(cdi = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi(t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "../cached_intermediates/1_acc_cdi_boot.rds")

accs_boot_cdi_byage <- acc_params |>
  partition(cluster) |>
  mutate(cdi = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi_age(t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(accs_boot_cdi_byage, "../cached_intermediates/1_acc_cdi_boot_byage.rds")
