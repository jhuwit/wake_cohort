# install.packages("fuzzyjoin")
# library(fuzzyjoin)
rm(list = ls()) 
# TAKE INDIVIDUAL FILES AND CREATE ONE DF WITH ALL TIME SERIES DATA 
check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Example usage:
required_packages <- c("readr", "dplyr", "tidyverse",
                       "lubridate", "purrr", "stringr", "here",
                       "slider", "furrr", "future", "patchwork", "zoo", "fuzzyjoin")
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
library(fuzzyjoin)
options(digits.secs = 2)
force = TRUE

ncores = parallelly::availableCores() - 1
options(dplyr.summarise.inform=FALSE)

covars = readxl::read_excel(
  here::here("data", "raw", "WAKE.flatfile.IDsfixed.11.27.24.xlsx")) %>% 
  mutate(ID = as.character(ID)) %>% 
  janitor::clean_names()

### to do: deal with ppl with >= 1 procedure 
covars %>% count(id) %>% filter(n!=1)

or_times = read_csv(here::here("data", "raw", "WAKE_OR_Times.IDfixed.11.29.24.csv")) %>% 
  janitor::clean_names()


or_times_small = 
  or_times %>% 
  select(id, anes_ready = anesthesia_ready, 
         anes_start = anesthesia_start,
         anes_end = anesthesia_stop,
         cpb_start = cv_bypass_initiated,
         cpb_end = cv_bypass_ended) %>% 
  mutate(across(-id, ~mdy_hm(.x, tz = "UTC")))

or_ids = unique(or_times_small$id)
covar_ids = unique(covars$id)

sum(!(covar_ids %in% or_ids))
# 30 pts in covar IDS not in or ids 
sum(!(or_ids %in% covar_ids))
# 552 pats in OR ids not in covar ids 

individual_files = list.files(here::here("data", "raw", "individual"),
                              recursive = TRUE,
                              full.names = TRUE)


ts_ids = sub(".*individual\\/(.+).csv.*", "\\1", individual_files)

sum(!(covar_ids %in% ts_ids))
sum(!(ts_ids %in% covar_ids))
sum(!(ts_ids %in% or_ids))



# read in all individual files 
if(!file.exists(here::here("data", "processed", "all_timeseries_stamped_raw.rds")) | force){
  # function to read in CI file 
  plan(multisession, workers = ncores)
  
  all_raw = 
    future_map_dfr(.x = individual_files,
                   .f = function(x){
                     read_csv(x) %>% 
                       select(timestamp, contains("anesthesia"), CI, MAP, CVP) %>% 
                       mutate(id = sub(".*individual\\/(.+).csv.*", "\\1", x))
                   }) %>% 
    janitor::clean_names()
  
  plan(sequential)
  if(!dir.exists(here::here("data", "processed"))){
    dir.create(here::here("data", "processed"), recursive = TRUE)
  }
  
  readr::write_rds(all_raw, 
                   here::here("data", "processed", "all_timeseries_stamped_raw.rds"), 
                   compress = "xz")
} 

hemo_ts_df = read_rds(here::here("data", "processed", "all_timeseries_stamped_raw.rds"))

# get time stamps 

### step two is to make full time series (from start of end to anesthesia) and then join with the hemodynamics data 
if(!file.exists(here::here("data", "processed", "hemo_timeseries.rds")) || force){
  ids = unique(hemo_ts_df$id)
  create_full_ts = 
    function(subject) {
      x = try({
        # get row for the subject
        times = or_times_small %>%
          filter(id == subject) 
        ts = hemo_ts_df %>% 
          filter(id == subject)
        # time sequence from start to end of anesthesia by minute
        timeseq =
          seq(times$anes_start,
              times$anes_end,
              by = "1 min")
        
        
        start_cpb = 
          times$cpb_start
        
        end_cpb =
          times$cpb_end
        
        
        df = tibble(id = rep(subject, length(timeseq)),
                    time = timeseq) %>%
          mutate(
            CPB = case_when(
              time < start_cpb ~ "pre",
              time > end_cpb ~ "post",
              time >= start_cpb &
              time <= end_cpb ~ "intra"
            )
          ) %>% 
          left_join(ts %>% select(id, timestamp, ci, map, cvp), by = c("time" = "timestamp", "id"))
        df
      })
      x
    }
  
  # create the "full" time series data set 
  plan(multisession, workers = ncores)
  full_timeseries =
    future_map(.x = ids,
        .f = create_full_ts) 
  
  plan(sequential)
  
  full_timeseries_df = 
    keep(full_timeseries, ~ inherits(.x, "tbl_df")) %>% 
    bind_rows()
  
  
  readr::write_rds(full_timeseries_df, here::here("data", "processed", "hemo_timeseries.rds"), compress = "xz")
}


