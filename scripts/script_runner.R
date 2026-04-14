source("../helper/common.R")

list.files() |>
  as_tibble() |>
  filter(str_detect(value, "^[1-3]_")) |>
  mutate(foo = walk(value, \(v){
    message("Running: ", v)

    env <- new.env(parent = globalenv())
    sys.source(basename(v), envir = env)
    rm(env)
    gc()
    message("Done: ", basename(v))
  }))
