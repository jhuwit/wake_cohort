
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
                       "slider", "furrr", "future", "patchwork", "zoo", "janitor")
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
rm(list = ls())
force = FALSE

source(here::here("code", "processing_fns.R"))
source(here::here("code", "filtering_settings.R"))

ncores = parallelly::availableCores() - 1

options(dplyr.summarise.inform=FALSE)

# read in covariates 
covars = readxl::read_excel(
  here::here("data", "raw", "WAKE.flatfile.IDsfixed.11.27.24.xlsx")) %>% 
  mutate(ID = as.character(ID)) %>% 
  janitor::clean_names()

or_times = read_csv(here::here("data", "raw", "WAKE_OR_Times.IDfixed.11.29.24.csv")) %>% 
  janitor::clean_names()


start_ends = 
  or_times %>% 
  select(id, anes_ready = anesthesia_ready, 
         anes_start = anesthesia_start,
         anes_end = anesthesia_stop,
         cpb_start = cv_bypass_initiated,
         cpb_end = cv_bypass_ended) %>% 
  mutate(across(-id, ~mdy_hm(.x, tz = "UTC")))


# read in time series 
full_timeseries_df = readRDS(here::here("data", "processed", "hemo_timeseries.rds"))


# add in "anes ready times" 
anes_ready = 
  start_ends %>% 
  select(id, time_anes_ready = anes_ready)

# now add in ANES ready times 
all_hemo = 
  full_timeseries_df %>% 
  left_join(anes_ready, by = "id") %>% 
  mutate(bin_anes_ready = if_else(time < time_anes_ready, 0, 1)) %>% 
  mutate(across(c(map, ci, cvp), ~if_else(is.nan(.x), NA_real_, .x))) %>% 
  select(-time_anes_ready)


# look at only post anes-ready measure and cut off 5 minutes before last non-NA map 
all_hemo = 
  all_hemo %>% 
  group_by(id) %>% 
  mutate(lastmap = max(time[!is.na(map)]) - as.period(5, "minutes"))

all_hemo_postready = 
  all_hemo %>% 
  filter(bin_anes_ready == 1 & time <= lastmap)


all_hemo_postready %>% 
  group_by(id) %>% 
  summarize(time = difftime(max(time), min(time), units = "hours")) %>% 
  arrange(time) %>% 
  tail()


all_hemo_postready %>% 
  group_by(id) %>% 
  summarize(time = difftime(max(time), min(time), units = "hours")) %>% 
  arrange(time)

length(unique(all_hemo_postready$id)) # how many ids we have 

summary(all_hemo_postready$ci, na.rm = TRUE)
nrow(all_hemo_postready %>% filter(ci <= 0.5)) / nrow(all_hemo_postready %>% filter(!is.na(ci)))
nrow(all_hemo_postready %>% filter(ci <= 1)) / nrow(all_hemo_postready %>% filter(!is.na(ci)))
nrow(all_hemo_postready %>% filter(ci <= 1.5)) / nrow(all_hemo_postready %>% filter(!is.na(ci)))
nrow(all_hemo_postready %>% filter(ci <= 2)) / nrow(all_hemo_postready %>% filter(!is.na(ci)))



sub_list = split(all_hemo_postready, f = all_hemo_postready$id)


if(!dir.exists(here::here("data", "processed", "interp_individual"))){
  dir.create(here::here("data", "processed", "interp_individual"), recursive = TRUE)
}

# subj = sub_list %>% bind_rows() %>% filter(id == "49184699858408")
for(subj in sub_list){
  id = subj$id[1]
  out_file = file.path(here::here("data", "processed", "interp_individual", paste0(id, ".rds")))
  if(!file.exists(out_file) || force){
    print(out_file)
    x = try({
      subj = 
        subj %>% 
        ungroup() 
      # sliding change windows 
      # abs_change_df = 
      #   sliding_abs_change(df = subj,
      #                      window_size = 1,  
      #                      time_var = "time",
      #                      hemo_vars = "cvp") # don't need a function for this since window is only 1 minute
      
      spike_df = 
        sliding_spike(df = subj,
                      window_size = 2, # whatever corresponds to two minutes 
                      time_var = "time",
                      hemo_vars = c("cvp"))
      # join sliding dfs with original data frame and ID flags 
      joined = 
        subj %>% 
        select(time, map, cvp, ci, CPB) %>% 
        # left_join(abs_change_df, by = "time") %>% 
        left_join(spike_df, by = "time") %>% 
        mutate(map_abs_change = abs(map - lag(map, n = 1L)),
               cvp_abs_change = abs(cvp - lag(cvp, n = 1L))) %>% 
        # now impose flags 
        mutate(
          cvp_flag = case_when(
            is.na(cvp) ~ NA_character_,
            cvp > cvp_max_thresh | cvp < cvp_min_thresh ~ "Magnitude flag",
            cvp_abs_change > cvp_abs_diff_thresh ~ "Absolute change flag",
            cvp_change > cvp_spike_thresh & cvp_post <= (cvp_baseline * cvp_spike_rtb_thresh) ~ "Spike flag",
            .default = "No flag"
          ),
          map_flag = case_when(
            is.na(map) ~ NA_character_,
            map < map_min_thresh | map > map_max_thresh ~ "Magnitude flag",
            # is.na(ci) & CPB != "intra" ~ "NA ci flag",
            map_abs_change > map_abs_diff_thresh ~ "Absolute change flag",
            .default = "No flag"),
          ci_flag = case_when(
            is.na(ci) ~ NA_character_,
            # ci < 1 ~ "Magnitude flag",
            .default = "No flag"
          ))
      
     
      
      # find areas to remove 
      remove_cvp = 
        try({
          joined %>% 
            filter(!is.na(cvp_flag) & cvp_flag != "No flag") %>% 
            group_by(time, cvp_flag) %>% 
            mutate(
              flag_start = case_when(
                cvp_flag == "Magnitude flag" ~ time - as.period(cvp_mag_rm_pre, "second"),
                cvp_flag == "Absolute change flag" ~ time - as.period(cvp_abs_rm_pre, "second"),
                cvp_flag == "Spike flag" ~ time - as.period(cvp_spike_rm_pre, "second"),
              ),
              flag_end = case_when(
                cvp_flag == "Magnitude flag" ~ time + as.period(cvp_mag_rm_post, "second"),
                cvp_flag == "Absolute change flag" ~ time + as.period(cvp_abs_rm_post, "second"),
                cvp_flag == "Spike flag" ~ time + as.period(cvp_spike_rm_post, "second")
              )) %>%
            select(time, cvp_flag, flag_start, flag_end) %>% 
            expand(ts = seq(flag_start, flag_end, by = "1 min")) %>% 
            ungroup() %>% 
            select(ts) %>% 
            distinct() %>% 
            mutate(cvp_flag = TRUE)
        })
      
      
      remove_map = 
        try({
          joined %>% 
            filter(!is.na(map_flag) & map_flag != "No flag") %>% 
            group_by(time, map_flag) %>% 
            mutate(
              flag_start = case_when(
                map_flag == "Magnitude flag" ~ time - as.period(map_mag_rm_pre, "second"),
                # map_flag == "NA ci flag" ~ time,
                map_flag == "Absolute change flag" ~ time - as.period(map_abs_rm_pre, "second")
              ),
              flag_end = case_when(
                map_flag == "Magnitude flag" ~ time + as.period(map_mag_rm_post, "second"),
                # map_flag == "NA ci flag" ~ time,
                map_flag == "Absolute change flag" ~ time + as.period(map_abs_rm_post, "second")
              )
            ) %>%
            expand(ts = seq(flag_start, flag_end, by = "1 min")) %>%
            ungroup() %>% 
            select(ts) %>% 
            distinct() %>% 
            mutate(map_flag = TRUE)
        })
      
      filtered = 
        subj %>% 
        left_join(if (is.data.frame(remove_cvp)) remove_cvp else tibble(ts = subj$time, cvp_flag = FALSE), by = c("time" = "ts")) %>%
        left_join(if (is.data.frame(remove_map)) remove_map else tibble(ts = subj$time, map_flag = FALSE), by = c("time" = "ts")) %>%
        mutate(
          cvp_filt =
            case_when(is.na(cvp) ~ NA_real_,
                      cvp_flag ~ NA_real_,
                      TRUE ~ cvp),
          map_filt =
            case_when(is.na(map) ~ NA_real_,
                      map_flag ~ NA_real_,
                      TRUE ~ map),
          ci_filt = if_else(ci < ci_min_thresh, 1, ci))
      
      
      ### next interpolation 
      map_interp =
        try({run_interpolation(filtered, hemo_var = "map_filt", time_var = "time", rle_thresh = interp_thresh)})
      
      cvp_interp =
        try({
        run_interpolation(filtered, hemo_var = "cvp_filt", time_var = "time", rle_thresh = interp_thresh)
        })
      
      ci_interp =
        try({run_interpolation(filtered, hemo_var = "ci_filt", time_var = "time", rle_thresh = interp_thresh)})
      
      hemo_interp =
        subj %>% 
        left_join(if (is.data.frame(map_interp)) map_interp else tibble(time = filtered$time, map_filt = filtered$map_filt, map_filt_interp = NA), by =c("time" = "time")) %>% 
        left_join(if (is.data.frame(cvp_interp)) cvp_interp else tibble(time = filtered$time, cvp_filt = filtered$cvp_filt, cvp_filt_interp = NA), by =c("time" = "time")) %>% 
        left_join(if (is.data.frame(ci_interp)) ci_interp else tibble(time = filtered$time, ci_filt = filtered$ci_filt, ci_filt_interp = NA), by =c("time" = "time")) %>% 
        # left_join(cvp_interp, by = c("time" = "time")) %>% 
        # left_join(ci_interp, by = c("time" = "time")) %>% 
        select(id, time, cat_cpb = CPB, 
               map_filt, cvp_filt, ci_filt,
               val_map = map_filt_interp, val_ci = ci_filt_interp, 
               val_cvp = cvp_filt_interp)
      write_rds(hemo_interp, out_file, compress = "xz")
    })
    rm(x)
  }
}


if(!file.exists(here::here("data", "processed", "hemo_timeseries_interp.rds")) || force){
  files = list.files(here::here("data", "processed", "interp_individual"),
                     full.names = TRUE,
                     recursive = TRUE)
  all_ts = map_dfr(.x = files,
                   .f = readRDS)
  write_rds(all_ts, here::here("data", "processed", "hemo_timeseries_interp.rds"),
          compress = "xz")
}


