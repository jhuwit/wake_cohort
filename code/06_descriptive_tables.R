
library(here)
library(ggpubr)
library(tidyverse)
library(tidymodels)
library(patchwork)
library(gt) 
library(gtsummary)
library(paletteer)
theme_set(theme_light())
force = FALSE
options(dplyr.summarise.inform=F)
col1 = "#F06719FF"; col2 = "#1BA3C6FF"
col1a = "#7873C0FF"; col2a = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
source(here::here("code", "utilities.R"))

covars = read_rds(here::here("data", "processed", "covars_proc.rds"))
map_ci = read_rds(here::here("data", "analysis", "mapci_ranges.rds"))


wake_covars = 
  covars %>% 
  filter(id %in% map_ci$id) %>% 
  rename(bin_aki48h = bin_aki) %>% 
  select(id, val_age, cat_gender, val_bmi, val_egfr,
         bin_htn, bin_diabetes, bin_mi,
         bin_chf, cat_ef,
         bin_iabp, bin_stroke, bin_pvd,
         bin_cld, bin_betablocker, bin_acearb,
         bin_statin, val_creatlst, bin_redo, bin_emergent,
         bin_aki48h, euro_predmort, val_proctime, val_crystalloid,
         cat_rbc, val_cpbtime,val_hematocrit) %>% 
  mutate(cohort = "WAKE") %>% 
  mutate(val_ef = 
         case_when(cat_ef == "20-25%" ~ 22.5,
                   cat_ef == "40-45%" ~ 42.5,
                   cat_ef == "45-50%" ~ 47.5,
                   cat_ef == "65-70%" ~ 67.5,
                   cat_ef == "15-20%" ~ 17.5,
                   cat_ef == "30-35%" ~ 32.5,
                   cat_ef == "35-40%" ~ 37.5,
                   cat_ef == "25-30%" ~ 27.5,
                   cat_ef == "70-75%" ~ 77.5,
                   cat_ef == "50-55%" ~ 52.5,
                   cat_ef == "60-65%" ~ 62.5,
                   cat_ef == "55-60%" ~ 57.5)) %>% 
  mutate(bin_ef40 = case_when(
    cat_ef %in% c("15-20%", "20-25%", "25-30%", "30-35%", "35-40%") ~ "<=40%",
    cat_ef %in% c("40-45%" ,"45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%") ~ ">40%",
    .default = NA_character_)) 

wake_covars =
  wake_covars %>% 
  drop_na()

map_ci = map_ci %>% 
  filter(id %in% wake_covars$id) 

hemo = read_rds(here::here("data", "processed", "hemo_analytic.rds"))  %>% 
  filter(!is.na(cat_cpb)) %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "intra", "post"))) %>% 
  filter(id %in% wake_covars$id) 


map_only = 
  hemo %>% 
  group_by(id) %>% 
  summarize(map_65 = sum(val_map < 65, na.rm = TRUE),
            .groups = "drop")

four_cats = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(low_map_low_ci = sum(val_map < 65 & val_ci < 2, na.rm = TRUE),
            low_map_hi_ci = sum(val_map < 65 & val_ci >= 2, na.rm = TRUE),
            hi_map_low_ci = sum(val_map >= 65 & val_ci < 2, na.rm = TRUE),
            hi_map_hi_ci = sum(val_map >= 65 & val_ci >= 2, na.rm = TRUE),
            .groups = "drop")



## map ci 


ci_ranges = read_rds(here::here("data", "analysis", "ci_ranges.rds")) %>% 
  filter(id %in% wake_covars$id)
ci_data = read_rds(here::here("data", "analysis", "ci_data.rds")) %>% 
  filter(id %in% wake_covars$id)


map_ci_post = read_rds(here::here("data", "analysis", "mapci_ranges_post.rds")) %>% 
  filter(id %in% wake_covars$id)

map_ci_pre = read_rds(here::here("data", "analysis", "mapci_ranges_pre.rds")) %>% 
  filter(id %in% wake_covars$id)


map_ci_cpb = 
  map_ci_post %>% mutate(cat_cpb = "post") %>% 
  bind_rows(map_ci_pre %>% mutate(cat_cpb = "pre")) 

ci_ranges4 = read_rds(here::here("data", "analysis", "ci_ranges4.rds")) %>% 
  filter(id %in% wake_covars$id)

ci_ranges7 = read_rds(here::here("data", "analysis", "ci_ranges7.rds")) %>% 
  filter(id %in% wake_covars$id)


map_data = read_rds(here::here("data", "analysis", "map_ranges.rds")) %>% 
  filter(id %in% wake_covars$id)
mapci_tert = read_rds(here::here("data", "analysis", "mapci_ranges_tertile.rds")) %>% 
  filter(id %in% wake_covars$id)

mapci_quint = read_rds(here::here("data", "analysis", "mapci_ranges_quintile.rds")) %>% 
  filter(id %in% wake_covars$id)


## some stats for paper 

four_cats_cpb  %>% 
  group_by(cat_cpb) %>% 
  summarize(across(low_map_low_ci,
                   .fns = list(median = ~median(.x),
                                q1 = ~quantile(.x, .25), 
                                q3 = ~quantile(.x, .75)),
                   .names = "{.col}_{.fn}"))
hemo %>% 
  filter(cat_cpb != "intra") %>% 
  summarize(val = median(val_ci, na.rm = TRUE),
            q25 = quantile(val_ci, 0.25, na.rm = TRUE),
            q75 = quantile(val_ci, 0.75, na.rm = TRUE))


hemo %>% 
  filter(cat_cpb != "intra") %>%
  group_by(id) %>% 
  summarize(low_ci = sum(val_ci <= 2, na.rm = TRUE),
            .groups = "drop") %>% 
  summarize(med = median(low_ci),
            iqr = IQR(low_ci),
            q1 = quantile(low_ci, .25),
            q3 = quantile(low_ci, .75))

hemo %>% 
  filter(cat_cpb != "intra") %>%
  group_by(id) %>% 
  summarize(low_ci = sum(val_ci <= 2, na.rm = TRUE),
            .groups = "drop") %>% 
  filter(low_ci >= 5) 


missing = hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id, cat_cpb) %>% 
  summarize(start = min(time),
            end = max(time),
            n = n(),
            across(c(val_map, val_cvp, val_ci),
                   ~sum(is.na(.x))),
            .groups = "drop") %>% 
  mutate(elapsed = as.numeric(difftime(end, start), units = "mins")) %>% 
  mutate(across(c(val_map, val_cvp, val_ci), ~ .x / n * 100, .names = "pct_{col}")) %>% 
  group_by(cat_cpb) %>% 
  summarize(across(c(elapsed, contains("val")),
                   .fns = list(mean = mean, sd = sd)))
# write.csv(missing, here::here("manuscript", "missing_data_summary.csv"))


# Table 1

## By AKI

t = tableone::CreateTableOne(vars = colnames(wake_covars %>% select(-id, -euro_predmort, -cohort)),
                             strata = "bin_aki48h",
                             factorVars =  colnames(wake_covars %>% select(starts_with("bin"))),
                             data = wake_covars)

t 
t = print(t)
if (!file.exists(here::here("manuscript", "table_one.csv")) || force){
  write.csv(t, here::here("manuscript", "table_one.csv"))
}


ci_tertiles = 
  ci_data %>% 
  summarize(across(ci_20, ~quantile(.x, c(1/3, 2/3)))) %>% 
  pull(ci_20)
ci_tertiles 
temp_ci = 
  ci_data %>% 
  mutate(tertile = case_when(
    ci_20 <= ci_tertiles[1] ~ "t1",
    ci_20 <= ci_tertiles[2] ~ "t2",
    .default = "t3"
  )) %>% 
  select(tertile, id) %>% 
  left_join(wake_covars, by = "id")

t = tableone::CreateTableOne(vars = colnames(wake_covars %>% select(-id, -euro_predmort, -cohort)),
                             strata = "tertile",
                             factorVars =  colnames(wake_covars %>% select(starts_with("bin"))),
                             data = temp_ci)
t 
t = print(t)
if(!file.exists(here::here("manuscript", "ci_tertile_table.csv")) || force) write.csv(t, here::here("manuscript", "ci_tertile_table.csv"))


### Time with MAP < 65

map_65 = 
  hemo %>% 
  group_by(id) %>% 
  summarize(map_65 = sum(val_map < 65, na.rm = TRUE),
            .groups = "drop") 

map_60 = 
  hemo %>% 
  group_by(id) %>% 
  summarize(map_60 = sum(val_map < 60, na.rm = TRUE),
            .groups = "drop") 
map_65_cpb = 
  hemo %>% 
  group_by(id, cat_cpb) %>% 
  summarize(map_65 = sum(val_map < 65, na.rm = TRUE),
            .groups = "drop") 
map_tertiles = 
  map_65 %>% 
  summarize(across(map_65, ~quantile(.x, c(1/3, 2/3), na.rm = TRUE))) %>% 
  pull(map_65)

map_tertiles 

temp_map = 
  map_65 %>% 
  mutate(tertile = case_when(
    map_65 <= map_tertiles[1] ~ "t1",
    map_65 <= map_tertiles[2] ~ "t2",
    .default = "t3"
  )) %>% 
  select(tertile, id) %>% 
  left_join(wake_covars, by = "id")

t = tableone::CreateTableOne(vars = colnames(wake_covars %>% select(-id, -euro_predmort, -cohort)),
                             strata = "tertile",
                             factorVars =  colnames(wake_covars %>% select(starts_with("bin"))),
                             data = temp_map)
t
t = print(t)
if(!file.exists(here::here("manuscript", "map_tertile_table.csv")) || force) write.csv(t, here::here("manuscript", "map_tertile_table.csv"))



map_65_ci2 = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(map65_ci2 = sum(val_map < 65 & val_ci < 2, na.rm = TRUE),
            .groups = "drop") 

map_tertiles = 
  map_65_ci2 %>% 
  summarize(across(map65_ci2, ~quantile(.x, c(1/3, 2/3)))) %>% 
  pull(map65_ci2)
map_tertiles 

temp_map = 
  map_65_ci2 %>% 
  mutate(tertile = case_when(
    map65_ci2 <= map_tertiles[1] ~ "t1",
    map65_ci2 <= map_tertiles[2] ~ "t2",
    .default = "t3"
  )) %>% 
  select(tertile, id) %>% 
  left_join(wake_covars, by = "id")

t = tableone::CreateTableOne(vars = colnames(wake_covars %>% select(-id, -euro_predmort, -cohort)),
                             strata = "tertile",
                             factorVars =  colnames(wake_covars %>% select(starts_with("bin"))),
                             data = temp_map)
t
t = print(t)
if(!file.exists(here::here("manuscript", "mapci_tertile_table.csv")) || force) write.csv(t, here::here("manuscript", "mapci_tertile_table.csv"))


map_ci %>% 
  summarize(across(contains("map"), list(mean = mean, sd = sd, n = ~ sum(.x >= 5, na.rm = TRUE), pctn = ~sum(.x >= 5, na.rm = TRUE) / n()))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(map = sub(".*map_([^_]*)_.*", "\\1", name),
         ci = sub(".*ci_([^_]*)_.*", "\\1", name),
         metric = sub(".*_(.*)", "\\1", name)) %>% 
  select(value, map, ci, metric) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(pctn = pctn*100,
         across(c(mean, sd, pctn), ~sprintf("%.2g", signif(.x, 2)))) %>% 
  mutate(mean_sd = paste0(mean, " (", sd, ")"),
         n_pct = paste0(n, " (", pctn, ")")) %>% 
  mutate(map = c(rep("<65", 4), rep(">=65", 4)),
         ci = rep(c("<=2", "(2,2.4]", "(2.4, 2.8]", ">2.8"), 2)) %>% 
  select(map, ci, mean_sd, n_pct)  %>% 
  magrittr::set_colnames(c("MAP", "CI", "Mean (SD)", "n (%)")) %>% 
  gt::gt()




