library(here)
# library(ggpubr)
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


##### ----- loading data --------
source(here::here("code", "prep_data.R"))

##### regressions 

## Univariate regressions 

df = 
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars) 

q = quantile(map_auc$map_auc,.99) %>% unlist()

df = 
  df %>% 
  left_join(map_only, by = "id") %>% 
  left_join(map_auc, by = "id") %>% 
  mutate(map_auc_win = if_else(map_auc > q, q, map_auc))

df2 = 
  df %>% 
  left_join(four_cats, by = "id") %>% 
  select(-map_65)

### correlations, check for multicollinearity 

cor_df = df %>% 
  select(contains("q"), contains("map_65"))

cor_mat = cor(cor_df)

ggcorrplot::ggcorrplot(cor_mat,
           type = "lower",
           lab = TRUE)

cor_df %>% 
  mutate(sum = hypo_q1 + hypo_q2 + hypo_q3 + hypo_q4) %>% 
  select(sum, map_65)

################# ----- running models ------- ############
result_uni = 
  map_dfr(.x = df %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("bin_aki48h", "~", x))
            model = glm(formula, data = df %>% mutate(across(contains("q"), ~.x / 5)), family = binomial)
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(2)
          })

if(!file.exists(here::here("manuscript", "univariate_regressions.csv")) || force) {
  result_uni %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "univariate_regressions.csv"))
  
}




## fully adjusted 
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

car::vif(model)

if(!file.exists(here::here("manuscript", "adjusted_regressions.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "adjusted_regressions.csv"))
} 

## fully adjusted model with MAP < 60
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_60,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)

summary(df$map_60)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   14.00   22.00   26.55   35.00  179.00 
if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map60.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map60.csv"))
} 

# fully adjusted model with MAP 55
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_55,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
# 
# summary(df$map_55)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    5.00    9.00   11.65   15.00   95.00 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map55.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map55.csv"))
} 

# fully adjusted model with MAP AUC

# 
# summary(df$map_auc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0   204.8   319.0   378.0   472.2  2791.0 


model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_auc,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ),
  family = binomial()
)
# 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_mapauc.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_mapauc.csv"))
} 

# exclude individuals in top 1% of procedure time 

proc_q = quantile(df$val_proctime, 0.99); proc_q
df_small = df %>% 
  filter(val_proctime <= proc_q)
df %>% nrow() - (df_small %>% nrow())


model = glm(
  bin_aki48h ~ .,
  data = df_small %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


if(!file.exists(here::here("manuscript", "revision", "adjusted_regressions_exclude_long_proc.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regressions_exclude_long_proc.csv"))
} 

## Nested models 

# Comparing MAP <65 model with all MAP,CI combos and MAP <65
# 
# **MAP < 65 univariate model** 




m1 = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_hematocrit,
    val_crystalloid,
    val_cpbtime,
    cat_rbc
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
m2 = 
  glm(
    bin_aki48h ~ .,
    data = df %>% select(
      contains("q"),
      map_65,
      bin_aki48h,
      cat_gender,
      val_creatlst,
      val_age,
      val_bmi,
      bin_htn,
      bin_diabetes,
      bin_stroke,
      bin_ef40,
      bin_mi,
      bin_chf,
      bin_pvd,
      bin_cld,
      bin_betablocker,
      bin_acearb,
      bin_statin,
      bin_redo,
      bin_emergent,
      val_hematocrit,
      val_crystalloid,
      val_cpbtime,
      cat_rbc
    ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
    family = binomial()
  )

m3 = 
  glm(
    bin_aki48h ~ .,
    data = df %>% select(
      hypo_q1,
      map_65,
      bin_aki48h,
      cat_gender,
      val_creatlst,
      val_age,
      val_bmi,
      bin_htn,
      bin_diabetes,
      bin_stroke,
      bin_ef40,
      bin_mi,
      bin_chf,
      bin_pvd,
      bin_cld,
      bin_betablocker,
      bin_acearb,
      bin_statin,
      bin_redo,
      bin_emergent,
      val_hematocrit,
      val_crystalloid,
      val_cpbtime,
      cat_rbc
    ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
    family = binomial()
  )

## check VIF

m2 %>% summary()
car::vif(m2)
m3 %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE)

if(!file.exists(here::here("manuscript", "map65_regression_adj.csv")) || force) {
  m1 %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "map65_regression_adj.csv"))
}

if(!file.exists(here::here("manuscript", "map65_regression_adj_bricks.csv")) || force) {
  m2 %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "map65_regression_adj_bricks.csv"))
}

## map 65 unadjusted
model = 
  glm(bin_aki48h ~ map_65,
      data = df %>% select(map_65, bin_aki48h) %>% mutate(across(contains("map"), ~ .x / 5)),
      family = binomial()
  )

if(!file.exists(here::here("manuscript", "map65_regression_unadj.csv")) || force) { 
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "map65_regression_unadj.csv"))
  
}





# **ANOVA comparing full and reduced model** 
  

anova(m1, m2, test = "Chisq")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1      1248     911.80                     
# 2      1240     900.83  8   10.964   0.2037


# **Comparing MAP <65 model with low MAP low CI model**
  

m1 = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    val_hematocrit,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
m2 = 
  glm(
    bin_aki48h ~ .,
    data = df %>% select(
      hypo_q1,
      map_65,
      bin_aki48h,
      cat_gender,
      val_creatlst,
      val_age,
      val_bmi,
      bin_htn,
      bin_diabetes,
      bin_stroke,
      bin_ef40,
      bin_mi,
      bin_chf,
      bin_pvd,
      bin_cld,
      bin_betablocker,
      bin_acearb,
      bin_statin,
      bin_redo,
      bin_emergent,
      val_crystalloid,
      val_cpbtime,
      cat_rbc,
      val_hematocrit
    ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
    family = binomial()
  )

# m2 %>% 
  # jtools::summ(exp = TRUE)
anova(m1, m2, test = "Chisq")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1      1248      911.8                       
# 2      1247      907.7  1   4.0946  0.04302 *

## Substitution 

# Sub 5 minutes MAP < 65, CI in quartile 1 (<2) to 5 minutes MAP < 65, CI in quartile 2 (2-2.4)


m1 = 
  glm(
    bin_aki48h ~ .,
    data = df %>% select(
      contains("q"),
      bin_aki48h,
      cat_gender,
      val_creatlst,
      val_age,
      val_bmi,
      bin_htn,
      bin_diabetes,
      bin_stroke,
      bin_ef40,
      bin_mi,
      bin_chf,
      bin_pvd,
      bin_cld,
      bin_betablocker,
      bin_acearb,
      bin_statin,
      bin_redo,
      bin_emergent,
      val_crystalloid,
      val_cpbtime,
      cat_rbc,
      val_hematocrit
    ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
    family = binomial()
  )

k  = matrix(c(0, -1, 1, rep(0, length(coef(m1)) - 3)), 1)
contr =  multcomp::glht(m1, linfct = k)
summary(contr)

## 95 % CI 
pt_est = -0.17617
se = 0.08846
exp(pt_est)
# [1] 0.8384754
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1] 0.7050041 0.9972155


# Sub 5 minutes MAP < 65, CI < 2 to 5 minutes MAP < 65, CI >= 2



m1 = glm(
  bin_aki48h ~ low_map_low_ci+  low_map_hi_ci + 
    hi_map_low_ci + hi_map_hi_ci + 
    cat_gender+
    val_creatlst+
    val_age+
    val_bmi+
    bin_htn+
    bin_diabetes+
    bin_stroke+
    bin_ef40+
    bin_mi+
    bin_chf+
    bin_pvd+
    bin_cld+
    bin_betablocker+
    bin_acearb+
    bin_statin+
    bin_redo+
    bin_emergent + 
    val_crystalloid + 
    val_cpbtime + 
    cat_rbc+
    val_hematocrit,
  data = df2 %>%  mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
k  <- matrix(c(0, -1, 1, rep(0, length(coef(m1)) - 3)), 1)
contr <- multcomp::glht(m1, linfct = k)
summary(contr)

## 95 % CI 
pt_est = -0.13957  
se = 0.06226
exp(pt_est)
# [1] 0.8697321
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1] 0.7698192 0.9826125]


# Sub 5 minutes MAP < 65, CI < 2 with 5 minutes MAP >= 65, CI < 2


k  <- matrix(c(0, -1, 0, 1, rep(0, length(coef(m1)) - 4)), 1)
contr <- multcomp::glht(m1, linfct = k)
summary(contr)


pt_est = -0.13164
se =0.05857 
exp(pt_est)
# [1] 0.8766565
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))

# [1] 0.7815804 0.9832982

##### ------  Regressions: subgroup analyses ----- ##### 
## Adjusted with all bricks in model 

### LVEF groups 

low_lvef = wake_covars %>% 
  filter(bin_ef40 == "<=40%") %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% low_lvef) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 3),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("LVEF <=40, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

low_lvef_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "LVEF <= 40")

mod1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "low_lvef")


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% low_lvef)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

hi_lvef_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "LVEF > 40")


p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 3),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("LVEF >40, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

mod2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "normal_lvef")

if(!file.exists(here::here("manuscript", "lvef_sensitivity.csv")) || force) { 
  mod1 %>%
    bind_rows(mod2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "lvef_sensitivity.csv"))
}


### Shock groups 


shock = wake_covars %>% 
  filter(bin_iabp == 1 | bin_emergent == 1) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% shock) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    # bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    # bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "shock")

shock_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "Shock")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("Shock, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% shock)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    # bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    # bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
noshock_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "No shock")

p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("No shock, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(type = "no shock")

if(!file.exists(here::here("manuscript", "shock_sensitivity.csv")) || force) {
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "shock_sensitivity.csv"))
}




low_egfr = wake_covars %>% 
  filter(val_egfr < 60) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% low_egfr) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

low_egfr_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "eGFR < 60")

m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "low egfr")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("eGFR < 60, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% low_egfr)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

hi_egfr_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "eGFR >= 60")

m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "normal egfr")
p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("eGFR >= 60, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

if(!file.exists(here::here("manuscript", "egfr_sensitivity.csv")) || force) { 
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "egfr_sensitivity.csv"))
}


## anemia groups 
anemic = wake_covars %>% 
  filter(bin_anemia_hct == 1) %>% 
  pull(id)


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% anemic) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

anemic_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "Anemic")

m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "anemic")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("Anemic, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% anemic)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

not_anemic_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "Not Anemic")

m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "not anemic")
p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("Not Anemic, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

if(!file.exists(here::here("manuscript", "revision", "anemia_sensitivity.csv")) || force) { 
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "revision", "anemia_sensitivity.csv"))
}


htn_hist = wake_covars %>% 
  filter(bin_htn == 1) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% htn_hist) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    # bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

htn_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "History of hypertension")

m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "htn")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("History of hypertension, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% htn_hist)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    # bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

not_htn_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "No history of hypertension")

m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "no htn")

p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("No history of hypertension, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

if(!file.exists(here::here("manuscript", "revision/htn_sensitivity.csv")) || force) { 
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "revision/htn_sensitivity.csv"))
}

# CHF
chf_hist = wake_covars %>% 
  filter(bin_chf == 1) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% htn_hist) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    # bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

chf_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "History of CHF")

m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "chf")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("History of CHF, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% chf_hist)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    # bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

not_chf_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "No history of CHF")

m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "no chf")

p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("No history of CHF, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

if(!file.exists(here::here("manuscript", "revision/chf_sensitivity.csv")) || force) { 
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "revision/chf_sensitivity.csv"))
}

# surgery phase

df_temp =  
  map_ci_pre %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "pre cpb")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("Pre-CPB"))

# model %>% 
#   broom::tidy(exponentiate = TRUE, conf.int = TRUE)
df_temp =  
  map_ci_post %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
# 
# model %>% 
#   broom::tidy(exponentiate = TRUE, conf.int = TRUE)
m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate(type = "post cpb")
p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = "Post-CPB")

if(!file.exists(here::here("manuscript", "cpbphase_sensitivity.csv")) || force) { 
  
  m1 %>%
    bind_rows(m2) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(here::here("manuscript", "cpbphase_sensitivity.csv"))
}



### Tertiles 


df_temp =  
  mapci_tert %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3",
                           "normo_q1", "normo_q2", "normo_q3")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

if(!file.exists(here::here("manuscript", "tertile_regression_adj.csv")) || force) { 
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    write_csv(here::here("manuscript", "tertile_regression_adj.csv"))
}

# check adition of MAP 

df_temp =  
  mapci_tert %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3",
                           "normo_q1", "normo_q2", "normo_q3")) %>% 
  left_join(wake_covars) %>% 
  left_join(map_only)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    map_65, 
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
  family = binomial()
)

model %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) 
### Quintiles 


df_temp =
  mapci_quint %>%
  filter(id %in% wake_covars$id) %>%
  magrittr::set_colnames(
    c(
      "id",
      "hypo_q1",
      "hypo_q2",
      "hypo_q3",
      "hypo_q4",
      "hypo_q5",
      "normo_q1",
      "normo_q2",
      "normo_q3",
      "normo_q4",
      "normo_q5"
    )
  ) %>%
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


if(!file.exists(here::here("manuscript", "quintile_regression_adj.csv")) || force) { 
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    write_csv(here::here("manuscript", "quintile_regression_adj.csv"))

}

## check with MAP 
df_temp =
  mapci_quint %>%
  filter(id %in% wake_covars$id) %>%
  magrittr::set_colnames(
    c(
      "id",
      "hypo_q1",
      "hypo_q2",
      "hypo_q3",
      "hypo_q4",
      "hypo_q5",
      "normo_q1",
      "normo_q2",
      "normo_q3",
      "normo_q4",
      "normo_q5"
    )
  ) %>%
  left_join(wake_covars) %>% left_join(map_only)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    map_65, 
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(c(contains("map"), contains("q")), ~ .x / 5)),
  family = binomial()
)

model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) 
# Investigating MAP associations 



m0 = 
  glm(bin_aki48h ~ .,
      data = df %>% 
        select(map_65, bin_aki48h) %>% 
        mutate(across(contains("map"), ~ .x / 5)),
      family = binomial())

m01 =
  glm(bin_aki48h ~ .,
      data = df %>%
        select(map_65, bin_aki48h, val_proctime) %>%
        mutate(across(contains("map"), ~ .x / 5)),
      family = binomial())

m1 = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_proctime,
    val_cpbtime
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)

m2 = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    bin_redo,
    bin_emergent,
    val_proctime,
    val_crystalloid,
    val_cpbtime,
    cat_rbc
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)




jtools::export_summs(m0, m01, m1, m2, exp = TRUE,
                     model.names = c("Uni", "MAP+Proctime", "MAP+Proctime+Demo", "MAP+Proctime+Demo+Intraop"),
                     error_format = "(p = {p.value})",
                     coefs = c("MAP<65" = "map_65",
                               "CPB time" = "val_cpbtime",
                               "Proc time" = "val_proctime",
                               "Crystalloid" = "val_crystalloid",
                               "RBC 0" = "cat_rbc0",
                               "RBC 1" = "cat_rbc1"))



###### probabilities 

quants = quantile(df$hypo_q1, c(seq(0, 1, .25)))

df_quantile = df %>% 
  mutate(quantile = cut(hypo_q1, breaks = c(0, 1, 6, 114), include.lowest = TRUE,
                        right = FALSE))
df_quantile %>% count(quantile)
df_quantile %>% 
  group_by(exposure_mins = quantile) %>% 
  summarize(n = n(),
            aki = sum(bin_aki48h),
            .groups = "drop") %>% 
  mutate(prob = aki / n)

df %>% 
  mutate(exposure = cut(hypo_q1, breaks= c(0, 5, 10, 15, 25, Inf),
         include.lowest = TRUE, right = FALSE)) %>% 
  group_by(exposure) %>% 
  summarize(n = n(),
            aki = sum(bin_aki48h),
            .groups = "drop") %>% 
  mutate(prob = aki / n)

g5 = 
  df %>% 
  filter(hypo_q1 >= 5)

g5 %>% 
  summarize(aki = sum(bin_aki48h),
            n = n()) %>% 
  mutate(prob = aki / n)


g5 = 
  df %>% 
  filter(hypo_q1 >= 5)

g5 %>% 
  summarize(aki = sum(bin_aki48h),
            n = n()) %>% 
  mutate(prob = aki / n)
