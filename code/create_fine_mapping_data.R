library(here)
library(tidyverse)
options(dplyr.summarise.inform=F)


# read in hemo data 
hemo_data = readRDS(here::here("data", "processed", "hemo_analytic.rds")) %>% 
  ungroup() 


# source regression functions 
source(here::here("code", "utilities.R"))
if(!dir.exists(here::here("data", "analysis"))){
  dir.create(here::here("data", "analysis"), recursive = TRUE)
}

### map cvp ranges, overall and stratified 
map_cvp = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra"),
              thresh1 = seq(45, 115, 10),
              var1 = "val_map",
              thresh2 = seq(0, 20, 2),
              var2 = "val_cvp")

write_rds(map_cvp, here::here("data", "analysis", "mapcvp_ranges.rds"))

map_cvp_cilow = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra", val_ci < 2),
              thresh1 = seq(45, 115, 10),
              var1 = "val_map",
              thresh2 = seq(0, 20, 2),
              var2 = "val_cvp")

write_rds(map_cvp_cilow, here::here("data", "analysis", "mapcvp_ranges_lowci.rds"))

map_cvp_cihi = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra", val_ci >= 2),
              thresh1 = seq(45, 115, 10),
              var1 = "val_map",
              thresh2 = seq(0, 20, 2),
              var2 = "val_cvp")

write_rds(map_cvp_cihi, here::here("data", "analysis", "mapcvp_ranges_hici.rds"))

### map ranges, overall and stratified 
thresholds = seq(45, 115, 5)
map_data = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                     thresholds = thresholds,
                     hemo_variable = "map") %>% 
  select(-contains("missing"))

write_rds(map_data, here::here("data", "analysis", "map_ranges.rds"))



map_data_cilow = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra", val_ci < 2),
                      thresholds = thresholds,
                      hemo_variable = "map") %>% 
  select(-contains("missing"))

write_rds(map_data_cilow, here::here("data", "analysis", "map_ranges_lowci.rds"))


map_data_cihi = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra", val_ci >= 2),
                            thresholds = thresholds,
                            hemo_variable = "map") %>% 
  select(-contains("missing"))

write_rds(map_data_cihi, here::here("data", "analysis", "map_ranges_hici.rds"))

### cvp ranges, overall and stratified 

thresholds = seq(0, 20, 2)
cvp_data = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra"),
                      thresholds = thresholds,
                      hemo_variable = "cvp") %>% 
  select(-contains("missing"))

write_rds(cvp_data, here::here("data", "analysis", "cvp_ranges.rds"))



cvp_data_cilow = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra", val_ci < 2),
                            thresholds = thresholds,
                            hemo_variable = "cvp") %>% 
  select(-contains("missing"))

write_rds(cvp_data_cilow, here::here("data", "analysis", "cvp_ranges_lowci.rds"))


cvp_data_cihi = get_ranges(hemo_data = hemo_data %>% filter(cat_cpb != "intra", val_ci >= 2),
                           thresholds = thresholds,
                           hemo_variable = "cvp") %>% 
  select(-contains("missing"))

write_rds(cvp_data_cihi, here::here("data", "analysis", "cvp_ranges_hici.rds"))



range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}

### zones 

shell_df = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
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

write_rds(shell_df, here::here("data", "analysis", "zones_five.rds"))

shell_df_lowci = 
  hemo_data %>% 
  filter(cat_cpb != "intra", val_ci < 2) %>% 
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

write_rds(shell_df_lowci, here::here("data", "analysis", "zones_five_lowci.rds"))

shell_df_hici = 
  hemo_data %>% 
  filter(cat_cpb != "intra", val_ci >= 2) %>% 
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

write_rds(shell_df_hici, here::here("data", "analysis", "zones_five_hici.rds"))


shell_df_3 = 
  hemo_data %>% 
  filter(cat_cpb != "intra") %>% 
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

write_rds(shell_df_3, here::here("data", "analysis", "zones_three.rds"))

shell_df_3_lowci = 
  hemo_data %>% 
  filter(cat_cpb != "intra", val_ci < 2) %>% 
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

write_rds(shell_df_3_lowci, here::here("data", "analysis", "zones_three_lowci.rds"))

shell_df_3_hici = 
  hemo_data %>% 
  filter(cat_cpb != "intra", val_ci >= 2) %>% 
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

write_rds(shell_df_3_hici, here::here("data", "analysis", "zones_three_hici.rds"))


shell_df_3_hici


### join all datasets -- overall CI and stratified 

all_data = 
  cvp_data %>% 
  full_join(map_data, by = "id") %>% 
  full_join(map_cvp, by = "id") %>% 
  full_join(shell_df %>% rename_with(.cols = starts_with("zone"), ~sub("zone", "shell", .x)), by = "id") %>% 
  full_join(shell_df_3, by = "id")

write_rds(all_data, here::here("data", "analysis", "all_ranges.rds"))


all_data = 
  cvp_data %>% 
  full_join(map_data, by = "id") %>% 
  full_join(map_cvp, by = "id") %>% 
  full_join(shell_df %>% rename_with(.cols = starts_with("zone"), ~sub("zone", "shell", .x)), by = "id") %>% 
  full_join(shell_df_3, by = "id")

write_rds(all_data, here::here("data", "analysis", "all_ranges.rds"))


all_data_low = 
  cvp_data_cilow %>% 
  full_join(map_data_cilow, by = "id") %>% 
  full_join(map_cvp_cilow, by = "id") %>% 
  full_join(shell_df_lowci %>% rename_with(.cols = starts_with("zone"), ~sub("zone", "shell", .x)), by = "id") %>% 
  full_join(shell_df_3_lowci, by = "id")

extra_ids = tibble(id = all_data$id[!(all_data$id %in% all_data_low$id)])

all_data_low = 
  all_data_low %>% bind_rows(extra_ids) %>% 
  mutate(across(.cols = everything(), ~replace_na(.x, 0)))

write_rds(all_data_low, here::here("data", "analysis", "all_ranges_lowci.rds"))

all_data_hi = 
  cvp_data_cihi %>% 
  full_join(map_data_cihi, by = "id") %>% 
  full_join(map_cvp_cihi, by = "id") %>% 
  full_join(shell_df_hici %>% rename_with(.cols = starts_with("zone"), ~sub("zone", "shell", .x)), by = "id") %>% 
  full_join(shell_df_3_hici, by = "id")

extra_ids = tibble(id = all_data$id[!(all_data$id %in% all_data_hi$id)])

all_data_hi = 
  all_data_hi %>% bind_rows(extra_ids) %>% 
  mutate(across(.cols = everything(), ~replace_na(.x, 0)))

write_rds(all_data_hi, here::here("data", "analysis", "all_ranges_hici.rds"))





### new physiology based zones 
map_cvp = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra"),
              thresh1 = c(0, 65, Inf),
              var1 = "val_map",
              thresh2 = c(0, 10, 20),
              var2 = "val_cvp",
              r = FALSE)

write_rds(map_cvp, here::here("data", "analysis", "pzones.rds"))

map_cvp_cilow = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra", val_ci < 2),
              thresh1 = c(0, 65, Inf),
              var1 = "val_map",
              thresh2 = c(0, 10, 20),
              var2 = "val_cvp", 
              r = FALSE)

write_rds(map_cvp_cilow, here::here("data", "analysis", "pzones_lowci.rds"))

map_cvp_cihi = 
  get_ranges2(hemo_data %>% filter(cat_cpb != "intra", val_ci >= 2),
              thresh1 = c(0, 65, Inf),
              var1 = "val_map",
              thresh2 = c(0, 10, 20),
              var2 = "val_cvp", 
              r = FALSE)

write_rds(map_cvp_cihi, here::here("data", "analysis", "pzones_hici.rds"))
