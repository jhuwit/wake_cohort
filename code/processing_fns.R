library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(here)
library(slider)
library(furrr)
library(future)
library(patchwork)
options(digits.secs = 2)


# function that includes NAs in RLE 
rle2 = function (x)  {
  stopifnot("'x' must be a vector of an atomic type" = is.atomic(x))
  
  n <- length(x)
  if (n == 0L) {
    return(structure(list(
      lengths = integer(), values = x)
    ), class = 'rle')
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Where does next value not equal current value?
  # i.e. y will be TRUE at the last index before a change
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- (x[-1L] != x[-n])
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Since NAs are never equal to anything, NAs in 'x' will lead to NAs in 'y'.
  # These current NAs in 'y' tell use nothing - Set all NAs in y to FALSE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y[is.na(y)] <- FALSE
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # When a value in x is NA and its successor is not (or vice versa) then that
  # should also count as a value change and location in 'y' should be set to TRUE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- y | xor(is.na(x[-1L]), is.na(x[-n]))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Any TRUE locations in 'y' are the end of a run, where the next value
  # changes to something different
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  i <- c(which(y), n)
  
  structure(list(
    lengths = diff(c(0L, i)),
    values  = x[i]
  ), class = 'rle')
}

# get absolute difference in min and max values from vector 
get_ad = function(vec){
  if(sum(!is.na(vec)) <= 1){
    return(NA)
  } else{
    abs(max(vec, na.rm = TRUE) - min(vec, na.rm = TRUE))
  }
}
# get abs diff of all columns in data frame 
abs_diff = function(df){
  df %>% 
    summarize(across(.cols = everything(),
                     .fns = list(abs_change = ~get_ad(.x))))
}

# example 
# 
# df = tibble(y = c(1:10, 20:10))
# abs_diff(df)
# df = tibble(y1 = seq(1:20),
#             y2 = seq(1, 40, by = 2))
# abs_diff(df)

# function to identify spike 
# get baseline value from vector (smallest value before max) 
get_baseline = function(vec){
  if(sum(!is.na(vec)) <= 1){
    return(NA)
  } else { 
    min(vec[1:which.max(vec)], na.rm = TRUE) 
  }
}
# get change from baseline of vector 
get_change = function(vec){
  if(sum(!is.na(vec)) <= 1){
    return(NA)
  } else { 
    max(vec, na.rm = TRUE) - min(vec[1:which.max(vec)], na.rm = TRUE)
  }
}
# get post-max value of vector 
get_post = function(vec){
  if(sum(!is.na(vec)) <= 1){
    return(NA)
  } else { 
    min(vec[which.max(vec):length(vec)], na.rm = TRUE)
  }
}
# get baseline, change from baseline, and post values from all columns in data frame 
spike = function(df){
  df %>% 
    summarize(across(.cols = everything(),
                     .fns = list(baseline = ~get_baseline(.x),
                                 change = ~get_change(.x),
                                 post = ~get_post(.x))))
                     # .fns = list(baseline = ~min(.x[1:which.max(.x)], na.rm = TRUE),
                     #             change = ~max(.x, na.rm = TRUE) - min(.x[1:which.max(.x)], na.rm = TRUE),
                     #             post = ~min(.x[which.max(.x):length(.x)], na.rm = TRUE))))
}


# example
# df = tibble(y = c(1:10, 20:10))
# spike(df)
# df = tibble(y1 = seq(1:20),
#             y2 = seq(1, 40, by = 2))
# spike(df)

# identify flag function
# function to apply absolute diff function across sliding windows in df 
sliding_abs_change = 
  function(df, window_size, time_var, hemo_vars){
    end_num = nrow(df) - (window_size -1)
    slide_dfr(.x = df,
              .before = 0,
              .after = window_size,
              # this is applied to each window 
              # get abs diff of each hemo var in addition to the start time of the window 
              .f = ~ abs_diff(.x %>% select(all_of(hemo_vars))) %>% 
                mutate(time = min(.x %>% pull(all_of(time_var))))) %>% 
      slice(-(end_num : nrow(df)))
}

# example 
# temp_df =
#   tibble(var1 = rnorm(100, mean = 10, sd = 10),
#          var2 = rnorm(100, 10, 1),
#          time = seq.POSIXt(from = as.POSIXct("2020-10-10 00:00:00"), to = as.POSIXct("2020-10-20 00:00:00"), length.out = 100))
# sliding_abs_change(temp_df, window_size = 5, time_var = "time", hemo_vars = c("var1", "var2"))

# function to apply spike function across sliding windows in df 
sliding_spike = 
  function(df, window_size, time_var, hemo_vars){
    end_num = nrow(df) - (window_size -1)
    slide_dfr(.x = df,
              .before = 0,
              .after = window_size,
              # this is applied to each window 
              # get abs diff of each hemo var in addition to the start time of the window 
              .f = ~ spike(.x %>% select(all_of(hemo_vars))) %>% 
                mutate(time = min(.x %>% pull(all_of(time_var))))) %>% 
      slice(-(end_num : nrow(df)))
  }
# example
# temp_df =
#   tibble(var1 = rnorm(100, mean = 10, sd = 10),
#          var2 = rnorm(100, 10, 1),
#          time = seq.POSIXt(from = as.POSIXct("2020-10-10 00:00:00"), to = as.POSIXct("2020-10-20 00:00:00"), length.out = 100))
# sliding_spike(temp_df, window_size = 5, time_var = "time", hemo_vars = c("var1", "var2"))


# interpolation function
run_interpolation = 
  function(df, hemo_var, time_var, rollmean_window = 3, rle_thresh){
    # filter to between first and last non-NA measurements 
    df_interp = 
      df %>% 
      select(time = all_of(time_var),
             hemo = all_of(hemo_var)) %>% 
      mutate(first_meas = min(time[!is.na(hemo)]),
             last_meas = max(time[!is.na(hemo)])) %>% 
      filter(between(time, first_meas, last_meas)) %>%
      select(-c(first_meas, last_meas)) %>% 
      mutate(
        time_num = seq(1, n()), # x's for interpolation 
        roll = case_when( # conditional rolling mean - previous 3 values if before NA, next 3 if after NA, NA if NA 
          is.na(hemo) ~ NA_real_, 
          is.na(lead(hemo, 1)) &
            !is.na(lag(hemo, 1)) ~ zoo::rollmean(
              hemo,
              k = rollmean_window,
              fill = NA,
              align = "right"
            ),
          !is.na(lead(hemo, 1)) &
            is.na(lag(hemo, 1)) ~ zoo::rollmean(
              hemo,
              k = rollmean_window,
              fill = NA,
              align = "left"
            )
        ),
        interp = approx( # do the interpolation 
          x = time_num,
          y = roll,
          xout = time_num,
          method = "linear"
        )$y
      )
    # now apply RLE to remove too long interpolations
    df_interp$rle = rep(rle2(df_interp$hemo)$lengths, times = rle2(df_interp$hemo)$lengths)
    
    # get rid of values if the length of missingness is > 15 min 
    df_interp %>% 
      mutate(final = if_else(is.na(hemo) & rle <= rle_thresh, interp, hemo)) %>% 
      select(time, hemo, final) %>% 
      magrittr::set_colnames(., c(time_var, hemo_var, paste0(hemo_var, "_interp")))
  }

                     
                        