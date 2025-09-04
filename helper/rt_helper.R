## RT function tested in peekbank_rts.Rmd
# takes rle_data dataframe (already rle'd)
get_rt <- function(rle_data, SAMPLING_RATE = 40, t_0 = TRUE, window_length = 0, t_end = FALSE, window_mostly_region = FALSE) {
  # end if no data
  print(rle_data)

  print(rle_data$values)

  if (is.null(rle_data$values) | is.null(rle_data$lengths)) {
    return(tibble(
      rt = NA,
      shift_type = NA
    ))
  }

  onset_aoi <- rle_data$values[1] # zero point AOI

  # end if missing for start
  if (t_0 & !(onset_aoi %in% c("target", "distractor"))) {
    return(tibble(
      rt = NA,
      shift_type = "onset_not_td"
    ))
  }


  first_departure_idx <- which(rle_data$values %in% c("target", "distractor"))[1]
  if (is.na(first_departure_idx)) {
    return(tibble(
      rt = NA,
      shift_type = "never_td"
    ))
  }

  onset_aoi <- rle_data$values[first_departure_idx] # not technically onset, but we're going with it for consistency

  other_aoi <- ifelse(onset_aoi == "distractor", "target", "distractor")

  if (window_length > 0 & (t_end | window_mostly_region)) {
    # TODO need to figure out from this data what is at end of window
    elapsed <- cumsum(rle_data$lengths)
    end_window_idx <- which(elapsed >= window_length)[1]
    if (is.na(end_window_idx)) {
      return(tibble(
        rt = NA,
        shift_type = "time_series_too_short"
      ))
    }
    if (other_aoi %in% rle$values[1:end_window_idx]) {
      return(tibble(
        rt = NA,
        shift_type = "both_td_in_window"
      ))
    }
    end_window_value <- rle_data$values[end_window_idx]
    if (t_end & end_window_value != onset_aoi) {
      return(tibble(
        rt = NA,
        shift_type = "window_end_not_td"
      ))
    }
    if (window_mostly_region) {
      in_region_idx <- which(rle_data$values[1:end_window_idx - 1] == onset_aoi) # indices in time region to relevant region
      if (length(in_region_idx) == 0) {
        time_points <- 0
      } else {
        time_points <- sum(rle_data$lengths[in_region_idx])
      }
      if (end_window_value == onset_aoi) { # count the partial amount
        time_points <- time_points + window_length - elapsed[end_window_idx - 1]
      }
      if (time_points < 3 / 4 * window_length) {
        return(tibble(
          rt = NA,
          shift_type = "window_not_enough_td"
        ))
      }
    }
  }


  first_landing <- rle_data$values[rle_data$values != onset_aoi &
    rle_data$values %in% c("target", "distractor")][1]

  # end if no shift
  if (is.na(first_landing)) {
    return(tibble(
      rt = NA,
      shift_type = "no shift"
    ))
  }

  shift_type <- case_when(
    onset_aoi == "distractor" &
      first_landing == "target" ~ "D-T",
    onset_aoi == "target" &
      first_landing == "distractor" ~ "T-D",
    TRUE ~ "other"
  )

  first_landing_idx <- which(rle_data$values == first_landing)[1]

  values_before_first_landing <- rle_data$values[1:(first_landing_idx - 1)]

  lengths_before_first_landing <- rle_data$lengths[1:(first_landing_idx - 1)]


  # we can either count from the first time they leave T/D or from the final time they leave T/D before landing on D/T
  # these are different in cases where they go D D D off off off D D T T ...

  last_shift_idx <- which(values_before_first_landing == onset_aoi)[length(which(values_before_first_landing == onset_aoi))]

  lengths_before_last_shift <- rle_data$lengths[1:last_shift_idx]

  length_before_first_shift <- rle_data$lengths[1]

  # if (first_landing_idx>2 & length(rle_data$lengths)>1) {
  #   shift_values <- rle_data$lengths[2:(first_landing_idx-1)]
  # } else {
  #   shift_values <- c()
  # }

  # rt is the number of samples happening before arrival + 1
  # (first sample of arrival)
  # times the length of a sample
  landing_time_rt <- (sum(lengths_before_first_landing) + 1) * (1000 / SAMPLING_RATE)

  # rt is the number of samples happening before a shift is initiated + 1
  # first sample off of onset AOI prior to moving to the next location
  # times the length of a sample
  shift_start_rt <- (length_before_first_shift + 1) * (1000 / SAMPLING_RATE)

  last_shift_rt <- (sum(lengths_before_last_shift) + 1) * (1000 / SAMPLING_RATE)

  # shift length: how long is the time between shift initiation and first landing
  shift_length <- landing_time_rt - shift_start_rt

  last_shift_length <- landing_time_rt - last_shift_rt

  return(tibble(
    rt = landing_time_rt, # define rt as landing time RT for now to avoid breaking things downstream
    shift_start_rt = shift_start_rt,
    shift_type = shift_type,
    shift_length = shift_length,
    last_shift_rt = last_shift_rt,
    last_shift_length = last_shift_length
  ))
}
