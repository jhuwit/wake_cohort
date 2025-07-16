library(here)
library(tidyverse)
options(dplyr.summarise.inform=F)


# read in hemo data 
hemo_data = readRDS(here::here("data", "processed", "hemo_analytic.rds")) %>% 
  ungroup() 

hemo_data %>% 
  group_by(cat_cpb, id) %>%
  summarize(ci_mins = sum(!is.na(val_ci))) %>% 
  group_by(cat_cpb) %>% 
  summarize(mean_mins = mean(ci_mins, na.rm = TRUE)) %>% 
  # filter(!is.na(val_cvp)) %>% 
  count() 

hemo_data %>% 
  filter(id == id[1]) %>% 
  pivot_longer(cols = c(val_map, val_ci, val_cvp)) 
# source regression functions 
source(here::here("code", "utilities.R"))
if(!dir.exists(here::here("data", "analysis"))){
  dir.create(here::here("data", "analysis"), recursive = TRUE)
}

#### data for analysis 

ci_vec = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
  pull(val_ci)

quantile(ci_vec, c(.25, .5, .75), na.rm = TRUE)


map_ci = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")
saveRDS(map_ci, here::here("data", "analysis", "mapci_ranges.rds"))

map_ci_pre = 
  get_ranges2(hemo_data %>% filter(cat_cpb == "pre"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")

saveRDS(map_ci_pre, here::here("data", "analysis", "mapci_ranges_pre.rds"))

map_ci_post = 
  get_ranges2(hemo_data %>% filter(cat_cpb == "post"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")
saveRDS(map_ci_post, here::here("data", "analysis", "mapci_ranges_post.rds"))

quantile(ci_vec, c(1/3, 2/3), na.rm = TRUE)

map_ci_tert = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2.2, 2.7, Inf),
              var2 = "val_ci")
saveRDS(map_ci_tert, here::here("data", "analysis", "mapci_ranges_tertile.rds"))

quantile(ci_vec, seq(0.2, 0.8, .2), na.rm = TRUE)

map_ci_quint = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 1.9, 2.2, 2.6, 2.9, Inf),
              var2 = "val_ci")
saveRDS(map_ci_quint, here::here("data", "analysis", "mapci_ranges_quintile.rds"))

ci_thresholds = seq(0, 10, 0.5)
ci_data = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                     thresholds = ci_thresholds,
                     hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data, here::here("data", "analysis", "ci_ranges.rds"))

ci_summary = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(ci_18 = sum(val_ci < 1.8, na.rm = TRUE),
            ci_16 = sum(val_ci< 1.6, na.rm = TRUE),
            ci_14 = sum(val_ci < 1.4, na.rm = TRUE),
            ci_12 = sum(val_ci < 1.2, na.rm = TRUE),
            ci_20 = sum(val_ci < 2, na.rm = TRUE),
            ci_22 = sum(val_ci < 2.2, na.rm = TRUE),
            ci_24 = sum(val_ci < 2.4, na.rm = TRUE),
            ci_26 = sum(val_ci < 2.6, na.rm = TRUE),
            ci_28 = sum(val_ci < 2.8, na.rm = TRUE),
            ci_30 = sum(val_ci < 3, na.rm = TRUE),
            minutes = n()) %>% 
  ungroup() %>% 
  mutate(across(starts_with("ci"), ~.x / minutes, .names = "{.col}_pct"))
saveRDS(ci_summary, here::here("data", "analysis", "ci_data.rds"))


map_thresholds = seq(45, 115, 5)
ci_thresholds = seq(0, 10, 0.5)

ci_thresh_4 = c(0, 2, 2.40, 2.80, Inf)
ci_thresh_7 = c(0, 1.8, 2.2, 2.40, 2.80, 3.20, 3.6, Inf)

# create time in range regression predictors 
map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "map") %>% 
  select(-contains("missing"))

saveRDS(map_data, here::here("data", "analysis", "map_ranges.rds"))



ci_data = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                     thresholds = ci_thresholds,
                     hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data, here::here("data", "analysis", "ci_ranges.rds"))

ci_data4 = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                      thresholds = ci_thresh_4,
                      hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data4, here::here("data", "analysis", "ci_ranges4.rds"))

ci_data7 = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                      thresholds = ci_thresh_7,
                      hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data7, here::here("data", "analysis", "ci_ranges7.rds"))



### 

#### 


#### old 
map_thresholds = seq(45, 115, 5)
map_thresholds_bivar = seq(45, 115, 10)
cvp_thresholds = seq(0, 20, 2)
ci_thresholds = seq(0, 10, 0.5)

ci_thresh_4 = c(0, 2, 2.40, 2.80, Inf)
ci_thresh_7 = c(0, 1.8, 2.2, 2.40, 2.80, 3.20, 3.6, Inf)

# create time in range regression predictors 
map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "map") %>% 
  select(-contains("missing"))

saveRDS(map_data, here::here("data", "analysis", "map_ranges.rds"))

cvp_data = get_ranges(hemo_data = hemo_data,
                      thresholds = cvp_thresholds,
                      hemo_variable = "cvp") %>% 
  select(-contains("missing"))

saveRDS(cvp_data, here::here("data", "analysis", "cvp_ranges.rds"))

ci_data = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                      thresholds = ci_thresholds,
                      hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data, here::here("data", "analysis", "ci_ranges.rds"))

ci_data4 = get_ranges(hemo_data = hemo_data,
                     thresholds = ci_thresh_4,
                     hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data4, here::here("data", "analysis", "ci_ranges4.rds"))

ci_data7 = get_ranges(hemo_data = hemo_data,
                      thresholds = ci_thresh_7,
                      hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data7, here::here("data", "analysis", "ci_ranges7.rds"))


ci_data4_cpb = get_ranges_cpb(hemo_data = hemo_data,
                      thresholds = ci_thresh_4,
                      hemo_variable = "ci") %>% 
  select(-contains("missing")) %>% 
  filter(cat_cpb != "intra")

saveRDS(ci_data4_cpb, here::here("data", "analysis", "ci_ranges4_cpb.rds"))

ci_data7 = get_ranges(hemo_data = hemo_data,
                      thresholds = ci_thresh_7,
                      hemo_variable = "ci") %>% 
  select(-contains("missing"))

saveRDS(ci_data7, here::here("data", "analysis", "ci_ranges7.rds"))

ci_data7_cpb = get_ranges_cpb(hemo_data = hemo_data,
                              thresholds = ci_thresh_7,
                              hemo_variable = "ci") %>% 
  select(-contains("missing")) %>% 
  filter(cat_cpb != "intra")

saveRDS(ci_data7_cpb, here::here("data", "analysis", "ci_ranges7_cpb.rds"))

bivar_bricks = get_ranges_biv(hemo_data, map_thresholds = map_thresholds_bivar,
                              cvp_thresholds = cvp_thresholds)
saveRDS(bivar_bricks, here::here("data", "analysis", "mapcvp_ranges.rds"))

# map-ci ranges

ci_vec = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
  pull(val_ci)

quantile(ci_vec, c(.25, .5, .75), na.rm = TRUE)


map_ci = 
  get_ranges2(hemo_data,
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")
saveRDS(map_ci, here::here("data", "analysis", "mapci_ranges.rds"))

map_ci_pre = 
  get_ranges2(hemo_data %>% filter(cat_cpb == "pre"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")

saveRDS(map_ci_pre, here::here("data", "analysis", "mapci_ranges_pre.rds"))

map_ci_post = 
  get_ranges2(hemo_data %>% filter(cat_cpb == "post"),
              thresh1 = c(0, 64, Inf),
              var1 = "val_map",
              thresh2 = c(0, 2, 2.4, 2.8, Inf),
              var2 = "val_ci")
saveRDS(map_ci_post, here::here("data", "analysis", "mapci_ranges_post.rds"))

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}

shell_df = 
  hemo_data %>% 
  mutate(
    shell = case_when(
      range_left(val_map, 95, 115) & between(val_cvp, 0, 8) ~ "zone_1",
      (range_left(val_map, 75, 95) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 8, 10)) ~ "zone_2",
      (range_left(val_map, 55, 75) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 65, 85) & range_left(val_cvp, 8, 12)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 10, 12)) ~ "zone_3",
      (between(val_map, 45, 55) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 55, 65) & range_left(val_cvp, 8, 12)) | 
        (range_left(val_map, 65, 115) & range_left(val_cvp, 12, 16)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 16, 18)) ~ "zone_4",
      (between(val_map, 45, 55) & range_left(val_cvp, 8, 20)) | 
        (range_left(val_map, 55, 65) & range_left(val_cvp, 12, 20)) | 
        (range_left(val_map, 65, 85) & range_left(val_cvp, 16, 20)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 18, 20)) ~ "zone_5"
    )
  ) %>% 
  group_by(id) %>% 
  count(shell, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell, values_from = n, id_cols = id) %>% 
  mutate(across(starts_with("zone"), ~ ifelse(is.na(.x), 0, .x)))

saveRDS(shell_df, here::here("data", "analysis", "zones_five.rds"))

shell_df_3 = 
  hemo_data %>% 
  mutate(
    shell = case_when(
      range_left(val_map, 95, 115) & between(val_cvp, 0, 8) ~ "group_1",
      (range_left(val_map, 75, 95) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 8, 10)) ~ "group_2",
      (range_left(val_map, 55, 75) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 65, 85) & range_left(val_cvp, 8, 12)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 10, 12)) ~ "group_3",
      (between(val_map, 45, 55) & between(val_cvp, 0, 8)) | 
        (range_left(val_map, 55, 65) & range_left(val_cvp, 8, 12)) | 
        (range_left(val_map, 65, 115) & range_left(val_cvp, 12, 16)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 16, 18)) ~ "group_4",
      (between(val_map, 45, 55) & range_left(val_cvp, 8, 20)) | 
        (range_left(val_map, 55, 65) & range_left(val_cvp, 12, 20)) | 
        (range_left(val_map, 65, 85) & range_left(val_cvp, 16, 20)) | 
        (range_left(val_map, 85, 115) & range_left(val_cvp, 18, 20)) ~ "group_5"
    )
  ) %>% 
  mutate(shell_new = 
           case_when(shell %in% c("group_1", "group_2") ~ "zone_1",
                     shell == "group_3" ~ "zone_2",
                     shell %in% c("group_4", "group_5") ~ "zone_3",
                     TRUE ~ NA)) %>% 
  select(-shell) %>% 
  group_by(id) %>% 
  count(shell_new, .drop = FALSE) %>% 
  ungroup() %>% 
  drop_na() %>% 
  pivot_wider(names_from =shell_new, values_from = n, id_cols = id) %>% 
  mutate(across(starts_with("zone"), ~ ifelse(is.na(.x), 0, .x)))

saveRDS(shell_df_3, here::here("data", "analysis", "zones_three.rds"))

map_thresholds = seq(45, 115, 10)
# map_thresholds_bivar = seq(45, 115, 10)
cvp_thresholds = seq(0, 20, 4)

# create time in range regression predictors 
map_data = get_ranges(hemo_data = hemo_data,
                      thresholds = map_thresholds,
                      hemo_variable = "map") %>% 
  select(-contains("missing"))

saveRDS(map_data, here::here("data", "analysis", "map_ranges_10.rds"))

cvp_data = get_ranges(hemo_data = hemo_data,
                      thresholds = cvp_thresholds,
                      hemo_variable = "cvp") %>% 
  select(-contains("missing"))

saveRDS(cvp_data, here::here("data", "analysis", "cvp_ranges_4.rds"))

### ci exposure 
ci_summary = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(ci_18 = sum(val_ci < 1.8, na.rm = TRUE),
            ci_16 = sum(val_ci< 1.6, na.rm = TRUE),
            ci_14 = sum(val_ci < 1.4, na.rm = TRUE),
            ci_12 = sum(val_ci < 1.2, na.rm = TRUE),
            ci_20 = sum(val_ci < 2, na.rm = TRUE),
            ci_22 = sum(val_ci < 2.2, na.rm = TRUE),
            ci_24 = sum(val_ci < 2.4, na.rm = TRUE),
            ci_26 = sum(val_ci < 2.6, na.rm = TRUE),
            ci_28 = sum(val_ci < 2.8, na.rm = TRUE),
            ci_30 = sum(val_ci < 3, na.rm = TRUE),
            minutes = n()) %>% 
  ungroup() %>% 
  mutate(across(starts_with("ci"), ~.x / minutes, .names = "{.col}_pct"))
saveRDS(ci_summary, here::here("data", "analysis", "ci_data.rds"))

ci_summary_strat = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id, cat_cpb) %>% 
  summarize(ci_18 = sum(val_ci < 1.8, na.rm = TRUE),
            ci_16 = sum(val_ci< 1.6, na.rm = TRUE),
            ci_14 = sum(val_ci < 1.4, na.rm = TRUE),
            ci_12 = sum(val_ci < 1.2, na.rm = TRUE),
            ci_20 = sum(val_ci < 2, na.rm = TRUE),
            ci_22 = sum(val_ci < 2.2, na.rm = TRUE),
            ci_24 = sum(val_ci < 2.4, na.rm = TRUE),
            ci_26 = sum(val_ci < 2.6, na.rm = TRUE),
            ci_28 = sum(val_ci < 2.8, na.rm = TRUE),
            ci_30 = sum(val_ci < 3, na.rm = TRUE),
            minutes = n()) %>% 
  ungroup() %>% 
  mutate(across(starts_with("ci"), ~.x / minutes, .names = "{.col}_pct"))
saveRDS(ci_summary_strat, here::here("data", "analysis", "ci_data_bycpb.rds"))


