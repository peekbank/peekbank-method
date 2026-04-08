source("../helper/common.R")

list.files() |>
  as_tibble() |>
  filter(str_detect(value, "^[1-4]_"), str_detect(value, "cdi|icc|test_retest")) |>
  mutate(foo = walk(value, \(v){
    message("Running: ", v)

    env <- new.env(parent = globalenv())
    sys.source(basename(v), envir = env)
    rm(env)
    gc()
    message("Done: ", basename(v))
  }))
