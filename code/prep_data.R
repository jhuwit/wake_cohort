

covars = read_rds(here::here("data", "processed", "covars_proc_extra.rds"))
map_ci = read_rds(here::here("data", "analysis", "mapci_ranges.rds"))

covars = covars %>% 
  filter(id %in% map_ci$id) %>% 
  mutate(bin_anemia_hgb = case_when(
    cat_gender == "Male" & val_hemoglobin < 13.5 ~ 1,
    cat_gender == "Male" & val_hemoglobin >= 13.5 ~ 0,
    cat_gender == "Female" & val_hemoglobin < 12 ~ 1,
    cat_gender == "Female" & val_hemoglobin >= 12 ~ 0,
  ),
  bin_anemia_hct = case_when(
    cat_gender == "Male" & val_hematocrit < 41 ~ 1,
    cat_gender == "Male" & val_hematocrit >= 41 ~ 0,
    cat_gender == "Female" & val_hematocrit < 36 ~ 1,
    cat_gender == "Female" & val_hematocrit >= 36 ~ 0,
  ))
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
         cat_rbc, val_cpbtime, val_hematocrit, bin_anemia_hct, bin_anemia_hgb,
         aki_stage, bin_aki7d) %>% 
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
            map_60 = sum(val_map < 60, na.rm = TRUE),
            map_55 = sum(val_map < 55, na.rm = TRUE),
            .groups = "drop")

four_cats = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(low_map_low_ci = sum(val_map < 65 & val_ci <= 2, na.rm = TRUE),
            low_map_hi_ci = sum(val_map < 65 & val_ci > 2, na.rm = TRUE),
            hi_map_low_ci = sum(val_map >= 65 & val_ci <= 2, na.rm = TRUE),
            hi_map_hi_ci = sum(val_map >= 65 & val_ci > 2, na.rm = TRUE),
            .groups = "drop")


four_cats_cpb = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id, cat_cpb) %>% 
  summarize(low_map_low_ci = sum(val_map < 65 & val_ci <= 2, na.rm = TRUE),
            low_map_hi_ci = sum(val_map < 65 & val_ci > 2, na.rm = TRUE),
            hi_map_low_ci = sum(val_map >= 65 & val_ci <= 2, na.rm = TRUE),
            hi_map_hi_ci = sum(val_map >= 65 & val_ci >2, na.rm = TRUE),
            .groups = "drop")

four_cats_cpb %>% 
  group_by(cat_cpb) %>% 
  summarize(med = median(low_map_low_ci),
            q1 = quantile(low_map_low_ci, 0.25),
            q3 =  quantile(low_map_low_ci, 0.75))

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

map_auc = read_rds(here::here("data", "analysis", "map_auc.rds"))
map_auc = 
  map_auc %>% 
  group_by(id) %>% 
  summarize(map_auc = sum(map_auc))
