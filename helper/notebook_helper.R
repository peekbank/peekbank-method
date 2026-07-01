library(tidyverse)
library(here)
library(tictoc)
library(ggthemes)
library(viridis)
library(ggh4x)
library(lme4)
library(broom.mixed)
library(cowplot)
library(emmeans)
knitr::opts_chunk$set(
  cache.extra = knitr::rand_seed, cache = FALSE,
  message = FALSE, warning = FALSE, error = FALSE, echo = F,
  fig.height = 4, fig.width = 6, dev = "cairo_pdf"
)
options(dplyr.summarise.inform = FALSE)

dataset_summ <- readRDS(here("cached_intermediates/paper_dataset_summary.rds"))

theme_set(theme_bw(base_size = 10))
theme_update(
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_line(),
  panel.grid.major.y = element_blank(),
  strip.background = element_blank(),
  legend.key = element_blank(),
  # panel.border = element_blank(),
  axis.line = element_line(),
  strip.text = element_text(face = "bold")
)

regroup_rt <- function(df) {
  df |>
    filter(time_0, time_end, during, frac == 1) |>
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
        str_detect(measure, "trim_last") ~ "trim_last",
        T ~ "untrimmed"
      )
    ) |>
    filter(trimming != "trim_last") |>
    ungroup() |>
    select(-time_0, -time_end, -during, -frac)
}


order_age <- function(df) {
  df |> mutate(
    age_bin = case_when(
      age_bin == "18-24" ~ "18–23 months",
      age_bin == "24-36" ~ "24–35 months",
      age_bin == "<18" ~ "<18 months",
      age_bin == ">=36" ~ ">=36 months"
    ),
    age_bin = factor(age_bin, levels = c("<18 months", "18–23 months", "24–35 months", ">=36 months"))
  )
}

reorder_within <- function(x, by, within) {
  new_x <- paste(x, within, sep = "___")
  stats::reorder(new_x, by)
}


scale_y_reordered <- function(...) {
  scale_y_discrete(labels = function(x) gsub("___.*$", "", x), ...)
}


filter_icc <- function(data, col_est) {
  include <- data |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    summarize(
      pct_0 = sum(.data[[col_est]] == 0) / n()
    )

  data |> semi_join(include |> filter(pct_0 < .5))
}

do_lollipop_plot <- function(df, name, include_condition, weights, flip = F, fixed_order = NULL) {
  df <- df |>
    mutate(name = name) |>
    {{ include_condition }}() |>
    filter(!is.na(combined)) |>
    left_join(weights) |>
    left_join(clean_names) |>
    group_by(`Dataset Name`) |>
    mutate(winner = if (flip) combined[which.min(est)] else combined[which.max(est)]) |>
    ungroup() |>
    mutate(`Dataset Name` = if (!is.null(fixed_order)) {
      factor(`Dataset Name`, levels = fixed_order)
    } else if (flip) {
      reorder(`Dataset Name`, -est)
    } else {
      reorder(`Dataset Name`, est)
    })

  p <- df |>
    ggplot(aes(x = `Dataset Name`, y = est)) +
    geom_line(aes(color = winner, group = `Dataset Name`), linewidth = 1.3, alpha = .5) +
    geom_point(aes(color = combined, size = n_admins)) +
    scale_size_area(max_size = 2) +
    guides(size = "none") +
    coord_flip() +
    labs(y = "Correlation") +
    scale_color_solarized(accent = "red") +
    facet_wrap(~name, nrow = 1) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(
      legend.position = "none", axis.title.y = element_blank(), # axis.text.y = element_text(size = 10)
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank()
    )

  if (!is.null(fixed_order)) p <- p + scale_x_discrete(drop = FALSE)
  p
}

do_lollipop_plot_age <- function(df, name, include_condition, weights, flip = F, xlim = NULL, fixed_order = NULL) {
  df <- df |>
    mutate(name = name) |>
    {{ include_condition }}() |>
    filter(!is.na(combined)) |>
    left_join(weights) |>
    left_join(clean_names) |>
    order_age() |>
    group_by(`Dataset Name`, age_bin) |>
    mutate(winner = if (flip) combined[which.min(est)] else combined[which.max(est)]) |>
    ungroup()

  if (!is.null(fixed_order)) {
    age_levels <- levels(df$age_bin)
    combined_levels <- c(outer(fixed_order, age_levels, paste, sep = "___"))

    # Add placeholder rows for missing dataset×age_bin combos so every facet has the same rows
    existing <- df |>
      distinct(`Dataset Name`, age_bin) |>
      filter(!is.na(`Dataset Name`))
    missing_combos <- crossing(`Dataset Name` = fixed_order, age_bin = age_levels) |>
      mutate(age_bin = factor(age_bin, levels = age_levels)) |>
      anti_join(existing, by = c("Dataset Name", "age_bin"))

    df <- bind_rows(df, missing_combos) |>
      mutate(
        age_bin = factor(age_bin, levels = age_levels),
        dataset_ordered = factor(paste(`Dataset Name`, age_bin, sep = "___"), levels = combined_levels)
      )
  } else {
    df <- df |>
      mutate(dataset_ordered = reorder_within(`Dataset Name`, if (flip) -est else est, age_bin))
  }

  missing_bins <- setdiff(levels(df$age_bin), as.character(unique(df$age_bin)))
  if (length(missing_bins) > 0) {
    df <- bind_rows(df, tibble(age_bin = factor(missing_bins, levels = levels(df$age_bin)), est = NA_real_))
  }

  df |>
    ggplot(aes(y = dataset_ordered, x = est)) +
    geom_line(aes(color = winner, group = `Dataset Name`), linewidth = 1.3, alpha = .5) +
    geom_point(aes(color = combined, size = n_admins)) +
    scale_size_area(max_size = 2) +
    guides(size = "none") +
    scale_y_reordered(expand = expansion(add = 1)) +
    coord_cartesian(xlim = xlim) +
    labs(title = name) +
    scale_color_solarized(accent = "red") +
    facet_wrap(~age_bin, scales = "free_y", ncol = 4, drop = FALSE, strip.position = "top") +
    annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf, lty = "dashed") +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      strip.placement = "outside",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8)
    )
}

set_of_age_lollis <- function(icc, include_func, is_rt = F, use_icc_order = FALSE) {
  fixed_order <- if (use_icc_order) {
    icc_means <- icc |>
      filter_icc("est") |>
      include_func() |>
      filter(!is.na(combined), !is.na(est)) |>
      left_join(clean_names) |>
      group_by(`Dataset Name`) |>
      summarize(mean_est = mean(est, na.rm = TRUE))

    icc |>
      distinct(dataset_name) |>
      left_join(clean_names) |>
      left_join(icc_means, by = "Dataset Name") |>
      mutate(mean_est = replace_na(mean_est, 0)) |>
      arrange(mean_est) |>
      pull(`Dataset Name`)
  } else {
    NULL
  }

  a <- do_lollipop_plot_age(
    icc |> filter_icc("est"), "ICC reliability", {{ include_func }},
    dataset_summ |> group_by(dataset_name) |> summarize(n_admins = n_distinct(administration_id)),
    fixed_order = fixed_order
  )
  g <- ggplotGrob(a + theme(legend.position = "bottom", legend.title = element_blank()))
  leg <- g$grobs[[which(g$layout$name == "guide-box-bottom")]]

  plot_grid(plot_grid(a, nrow = 1), leg, nrow = 2, rel_heights = c(1, .2))
}

set_of_lollis <- function(icc, cdi, trt, include_func, is_rt = F, use_icc_order = FALSE) {
  fixed_order <- if (use_icc_order) {
    icc_means <- icc |>
      filter_icc("est") |>
      include_func() |>
      filter(!is.na(combined), !is.na(est)) |>
      left_join(clean_names) |>
      group_by(`Dataset Name`) |>
      summarize(mean_est = mean(est, na.rm = TRUE))

    bind_rows(
      icc |> select(dataset_name),
      cdi |> select(dataset_name),
      trt |> select(dataset_name)
    ) |>
      distinct(dataset_name) |>
      left_join(clean_names) |>
      left_join(icc_means, by = "Dataset Name") |>
      mutate(mean_est = replace_na(mean_est, 0)) |>
      arrange(mean_est) |>
      pull(`Dataset Name`)
  } else {
    NULL
  }

  a <- do_lollipop_plot(
    icc |> filter_icc("est"), "ICC reliability", {{ include_func }},
    dataset_summ |> group_by(dataset_name) |> summarize(n_admins = n_distinct(administration_id)),
    fixed_order = fixed_order
  )
  b <- do_lollipop_plot(
    cdi |> filter(!is.na(comp_est)) |> rename(est = comp_est), "CDI comprehension", {{ include_func }},
    dataset_summ |> group_by(dataset_name) |> filter(!is.na(comp)) |> summarize(n_admins = n_distinct(administration_id)),
    flip = is_rt, fixed_order = fixed_order
  )
  c <- do_lollipop_plot(cdi |> filter(!is.na(prod_est)) |> rename(est = prod_est), "CDI production", {{ include_func }},
    dataset_summ |> group_by(dataset_name) |> filter(!is.na(prod)) |> summarize(n_admins = n_distinct(administration_id)),
    flip = is_rt, fixed_order = fixed_order
  )
  d <- do_lollipop_plot(
    trt |> filter(!is.na(est)), "Test-retest reliability", {{ include_func }},
    dataset_summ |> make_test_retest_pairs() |> group_by(dataset_name) |>
      summarize(n_admins = n_distinct(pair_number)),
    fixed_order = fixed_order
  )

  g <- ggplotGrob(a + theme(legend.position = "bottom", legend.title = element_blank()))
  leg <- g$grobs[[which(g$layout$name == "guide-box-bottom")]]

  inner <- plot_grid(plot_grid(a, d), NULL, plot_grid(b, c),
                     nrow = 1, rel_widths = c(1, 0.05, 1))
  labeled <- ggdraw(inner) +
    draw_label("A", x = 0,          y = 1, hjust = 0, vjust = 1.5, size = 14, fontface = "bold") +
    draw_label("B", x = 1.05 / 2.05, y = 1, hjust = 0, vjust = 1.5, size = 14, fontface = "bold")
  plot_grid(labeled, leg, nrow = 2, rel_heights = c(1, .1))
}

do_model <- function(df, make_baseline, weight_df, form) {
  data <- df |>
    filter(!is.na(est)) |>
    left_join(weight_df) |>
    mutate(across(-c(est, n_admins), as.factor)) |>
    {{ make_baseline }}()

  mod <- lmer(str_c("est ~", form, "+ (1 | dataset_name)"), data = data, weights = n_admins)
  emm <- emmeans(mod, as.formula(str_c("~", form)))
  contrast(emm, method = "trt.vs.ctrl", ref = 1) |> tidy(conf.int = T)
}


make_model_plot <- function(df, x, col = NULL, facet = NULL, fix_function, lab, breaks = c(-.1, 0, .1), limits = c(-.2, .2), dodge = 0, ref = NULL, ref_facet = NULL) {
  x_quo <- rlang::enquo(x)
  facet_quo <- rlang::enquo(facet)
  col_quo <- rlang::enquo(col)

  d <- df |>
    {{ fix_function }}()

  p <- d |>
    ggplot(aes(x = {{ x }}, col = {{ col }}, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange(position = position_dodge2(width = dodge)) +
    coord_flip() +
    labs(y = lab) +
    scale_y_continuous(breaks = breaks, limits = limits) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.position = "none", legend.title = element_blank())

  if (!rlang::quo_is_null(facet_quo)) {
    p <- p + facet_wrap(vars(!!facet_quo))
  }

  if (!rlang::quo_is_null(col_quo)) {
    ex <- rlang::get_expr(col_quo)
    col_name <- if (rlang::is_formula(ex)) {
      rlang::as_name(rlang::f_rhs(ex))
    } else {
      rlang::as_name(ex)
    }
    v <- d[[col_name]]

    if (is.factor(v) || is.character(v) || is.logical(v)) {
      p <- p + scale_color_solarized()
    } else if (is.numeric(v)) {
      p <- p + scale_color_viridis()
    } else {
      stop("Unsupported color column type")
    }
  }

  if (!is.null(ref)) {
    x_name <- rlang::as_name(x_quo)
    ref_df <- tibble(!!x_name := ref, estimate = 0)
    for (facet_col in names(ref_facet)) {
      ref_df[[facet_col]] <- ref_facet[[facet_col]]
    }
    p <- p + geom_point(
      data = ref_df, aes(x = .data[[x_name]], y = estimate),
      inherit.aes = FALSE, shape = 1, size = 3
    )
  }

  p
}

make_model_grid_plot <- function(icc, cdi, trt, make_baseline, form, x, col = NULL, facet = NULL, fix_function, breaks = c(-.1, 0, .1), limits = c(-.2, .2), dodge = 0, legend_height = .1, ref = NULL, ref_facet = NULL) {
  icc_plot <- do_model(
    icc |> filter_icc("est"), {{ make_baseline }},
    dataset_summ |> group_by(dataset_name) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |> make_model_plot({{ x }},
    col = {{ col }}, facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "ICC reliability", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )
  comp_plot <- do_model(
    cdi |> filter(!is.na(comp_est)) |> rename(est = comp_est), {{ make_baseline }},
    dataset_summ |> group_by(dataset_name) |> filter(!is.na(comp)) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |> make_model_plot({{ x }},
    col = {{ col }}, facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "Corr. w/ CDI Comprehension", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )
  prod_plot <- do_model(
    cdi |> filter(!is.na(prod_est)) |> rename(est = prod_est), {{ make_baseline }},
    dataset_summ |> group_by(dataset_name) |> filter(!is.na(prod)) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |> make_model_plot({{ x }},
    col = {{ col }}, facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "Corr. w/ CDI Production", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )
  trt_plot <- do_model(
    trt |> filter(!is.na(est)), {{ make_baseline }},
    dataset_summ |> make_test_retest_pairs() |> group_by(dataset_name) |>
      summarize(n_admins = n_distinct(pair_number)),
    form
  ) |> make_model_plot({{ x }},
    col = {{ col }}, facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "Test-retest reliability", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )

  g <- ggplotGrob(icc_plot + theme(legend.position = "bottom", legend.title = element_blank()))
  leg <- g$grobs[[which(g$layout$name == "guide-box-bottom")]]

  plot_grid(plot_grid(icc_plot, trt_plot), plot_grid(comp_plot, prod_plot), leg, nrow = 3, rel_heights = c(1, 1, legend_height))
}


bin_ages <- function(d_aoi, min_per_bin = 5) {
  d_aoi |>
    distinct(administration_id, age, dataset_name, across(any_of(c("prod", "comp", "pair_number", "est")))) |>
    mutate(age_bin = case_when(
      age < 18 ~ "<18",
      age < 24 ~ "18-24",
      age < 36 ~ "24-36",
      age >= 36 ~ ">=36"
    )) |>
    group_by(dataset_name, age_bin) |>
    mutate(count = n()) |>
    filter(count >= min_per_bin) |>
    ungroup()
}

safe_model_age <- function(d, form) {
  tryCatch(
    {
      mod <- lmer(str_c("est ~", form, "+ (1 | dataset_name)"), data = d, weights = n_admins)
      emm <- emmeans(mod, as.formula(str_c("~", form)))
      contrast(emm, method = "trt.vs.ctrl", ref = 1) |> tidy(conf.int = TRUE)
    },
    error = function(e) {
      message(
        "Failed for age_bin with ", nrow(d), " rows, ",
        n_distinct(d$dataset_name), " datasets: ", e$message
      )
      NULL
    }
  )
}
do_model_age <- function(df, make_baseline, weight_df, form) {
  data <- df |>
    filter(!is.na(est)) |>
    left_join(weight_df) |>
    mutate(across(!est & !n_admins, as.factor)) |>
    {{ make_baseline }}()

  data |>
    group_by(age_bin) |>
    nest() |>
    mutate(result = map(data, \(d){
      safe_model_age(d, form)
    })) |>
    select(-data) |>
    filter(!map_lgl(result, is.null)) |>
    unnest(result)
}


make_model_plot_age <- function(df, x, col = NULL, facet = NULL, fix_function, lab, breaks = c(-.1, 0, .1), limits = c(-.2, .2), dodge = .5, ref = NULL, ref_facet = NULL) {
  x_quo <- rlang::enquo(x)
  facet_quo <- rlang::enquo(facet)
  p <- df |>
    {{ fix_function }}() |>
    order_age() |>
    rename(Ages = age_bin) |>
    ggplot(aes(x = {{ x }}, col = Ages, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange(position = position_dodge2(width = dodge)) +
    coord_flip() +
    labs(y = lab) +
    scale_y_continuous(breaks = breaks, limits = limits) +
    geom_hline(yintercept = 0, lty = "dashed") +
    theme(legend.position = "none") +
    scale_color_manual(values = c(
      "<18 months" = "#FDE725", "18–23 months" = "#35B779",
      "24–35 months" = "#31688E", ">=36 months" = "#440154"
    ))

  if (!rlang::quo_is_null(facet_quo)) {
    p <- p + facet_wrap(vars(!!facet_quo))
  }

  if (!is.null(ref)) {
    x_name <- rlang::as_name(x_quo)
    ref_df <- tibble(!!x_name := ref, estimate = 0)
    for (facet_col in names(ref_facet)) {
      ref_df[[facet_col]] <- ref_facet[[facet_col]]
    }
    p <- p + geom_point(
      data = ref_df, aes(x = .data[[x_name]], y = estimate),
      inherit.aes = FALSE, shape = 1, size = 3
    )
  }

  p
}

make_model_grid_plot_age <- function(icc, cdi, make_baseline, form, x, facet = NULL, fix_function, dodge = .5, breaks = c(-.1, 0, .1), limits = c(-.2, .2), ref = NULL, ref_facet = NULL) {
  icc_plot <- do_model_age(
    icc |> filter_icc("est"), {{ make_baseline }},
    dataset_summ |> bin_ages() |> group_by(dataset_name, age_bin) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |>
    make_model_plot_age({{ x }},
      facet = {{ facet }}, fix_function = {{ fix_function }},
      lab = "ICC reliability", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
    )

  comp_plot <- do_model_age(
    cdi |> filter(!is.na(comp_est)) |> rename(est = comp_est), {{ make_baseline }},
    dataset_summ |> bin_ages() |> group_by(dataset_name, age_bin) |> filter(!is.na(comp)) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |> make_model_plot_age({{ x }},
    facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "Corr. w/ CDI Comprehension", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )
  prod_plot <- do_model_age(
    cdi |> filter(!is.na(prod_est)) |> rename(est = prod_est), {{ make_baseline }},
    dataset_summ |> bin_ages() |> group_by(dataset_name, age_bin) |> filter(!is.na(prod)) |> summarize(n_admins = n_distinct(administration_id)),
    form
  ) |> make_model_plot_age({{ x }},
    facet = {{ facet }}, fix_function = {{ fix_function }},
    lab = "Corr. w/ CDI Production", breaks = breaks, limits = limits, dodge = dodge, ref = ref, ref_facet = ref_facet
  )

  g <- ggplotGrob(icc_plot + theme(legend.position = "right"))
  leg <- g$grobs[[which(g$layout$name == "guide-box-right")]]

  plot_grid(plot_grid(icc_plot, leg), plot_grid(comp_plot, prod_plot), nrow = 2, rel_heights = c(1, 1))
}

clean_names <- tribble(
  ~dataset_name, ~`Dataset Name`, ~`Language`,
  "reflook_v4", "Yurovsky et al. (2017)", "English",
  "yurovsky_2017", "Yurovsky & Frank (2017)", "English",
  "weaver_zettersten_2024", "Weaver et al. (2024)", "English",
  "fernald_marchman_2012", "Fernald & Marchman (2012)", "English",
  "frank_tablet_2016", "Frank et al. (2016)", "English",
  "fmw_2013", "Fernald et al. (2013)", "English",
  "adams_marchman_2018", "Adams et al. (2018)", "English",
  "potter_canine", "Potter & Lew-Williams (2023)", "English",
  "fernald_totlot", "Fernald et al. (2006)", "English",
  "pomper_prime", "Pomper & Saffran (2015)", "English",
  "pomper_saffran_2016", "Pomper & Saffran (2016)", "English",
  "swingley_aslin_2002", "Swingley & Aslin (2002)", "English",
  "pomper_salientme", "Pomper & Saffran (2019)", "English",
  "perry_cowpig", "Perry & Saffran (2017)", "English",
  "ronfard_2021", "Ronfard et al. (2022)", "English",
  "bacon_gendercues", "Bacon & Saffran (2022)", "English",
  "garrison_bergelson_2020", "Garrison et al. (2020)", "English",
  "pomper_yumme", "Pomper & Saffran (2018)", "English",
  "mahr_coartic", "Mahr et al. (2015)", "English",
  "potter_remix", "Potter et al. (2019)", "English",
  "reflook_socword", "Yurovsky et al. (2013)", "English",
  "yoon_simpimp_2015", "Yoon et al. (2015)", "English",
  "borovsky_2019", "Borovsky & Peters (2019)", "English",
  "pomper_dimy", "Pomper & Saffran (2017)", "English",
  "weaver_saffran_2026", "Weaver & Saffran (2026)", "English",
  "bergelson_swingley_2012", "Bergelson & Swingley (2012)", "English",
  "moore_bergelson_2022_verb", "Moore & Bergelson (2022)", "English",
  "hurtado_2008", "Hurtado et al. (2008)", "Spanish",
  "kartushina_2019", "Kartushina & Mayor (2019)", "Norwegian",
  "weisleder_stl", "Weisleder & Fernald (2013)", "Spanish",
  "sander-montant_2022", "Sander-Montant et al. (2023)", "English, French",
  "gazetriggered_2020", "Egger et al. (2020)", "Dutch",
  "xsectional_2007", "Hurtado et al. (2007)", "Spanish",
  "casillas_tseltal_2015", "Casillas et al. (2017)", "Tseltal",
  "kremin_2021", "Kremin et al. (2023)", "English, French, Spanish",
  "byers-heinlein_2017", "Byers-Heinlein et al. (2017)", "English, French",
)


in_text_stats <- function(df, sign_flip = F, places = 2) {
  est <- df$estimate[1]
  low <- df$conf.low[1]
  high <- df$conf.high[1]
  if (sign_flip) {
    est <- -est
    new_low <- -high
    new_high <- -low
    low <- new_low
    high <- new_high
  }
  accuracy <- (1 / 10)**places
  str_c(scales::comma(est, accuracy), " [", scales::comma(low, accuracy = accuracy), ", ", scales::comma(high, accuracy), "]")
}
