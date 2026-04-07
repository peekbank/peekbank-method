library(tidyverse)
library(here)

list.files(here("scripts")) |> as_tibble() |>  filter(str_detect(value,"cdi|icc|test_retest")) |> mutate(foo=walk(value, \(v){
  message("Running: ",  v)
  
  env <- new.env(parent = globalenv())
  withr::with_dir(dirname("scripts"), {
    sys.source(basename(v), envir = env)
  })
  rm(env)
  gc()
  message("Done: ", basename(v))
}))
