
library(readr)
library(tidyverse)
library(purrr)
library(stringr)
library(lubridate)
library(dplyr)
library(tidyr)
library(mgcv)
library(here)
library(tidymodels)
`%notin%` <- Negate(`%in%`)
set.seed(123)
options(dplyr.summarise.inform=F)

normalize_quantiles = function(ref_data, new_data, col_name = "val_ci"){
  ref_vec = ref_data %>% pull(all_of(col_name)) %>% sort # automatically removes NAs 
  # ref_quantiles = seq(0, 1, length.out = length(ref_vec))
  key =  
    new_data %>%
    ungroup() %>% 
    mutate(row = row_number())
  
  df_sort = 
    new_data %>%
    ungroup() %>% 
    select(id, time, v = all_of(col_name)) %>% 
    mutate(row = row_number()) %>% 
    drop_na() %>% 
    arrange(v) 
  
  ci_vec = df_sort %>% pull(v)
  
  ci_vec_new = quantile(ref_vec, probs = ecdf(ci_vec)(ci_vec), type = 1) %>% unname
  # quantile(ci_vec, probs = seq(0,1, .1))
  # quantile(ci_vec_new, probs = seq(0,1, .1))
  # quantile(ref_vec, probs = seq(0,1, .1))
  
  res_tibble =
    tibble(new = ci_vec_new,
           row = df_sort$row)
  
  ret = 
    key %>% 
    left_join(res_tibble, by = "row") %>% 
    select(-row) %>% 
    rename(!!paste0(col_name, "_qnorm") := new)
  ret 
  
}

calc_euroscore = function(df){
  vars_nec = c("val_age", "cat_gender", "cr_g2_26", "bin_pvd", 
               "bin_cld", "bin_redo", "bin_mi", "ef_30_50", "ef_l30", 
               "bin_emergent")
  stopifnot(vars_nec %in% colnames(df))
  df %>% 
    mutate(gender = if_else(cat_gender == "Female", 1, 0)) %>% 
    mutate(z = (((1 + pmax(0, val_age - 59))* 0.0666354) + 
             (gender * 0.3304052) + 
             (cr_g2_26 * 0.6521653) + 
              (bin_pvd * 0.6558917) + 
    (bin_cld * 0.4931341)+ 
    (bin_redo * 1.002625) + 
    (bin_mi * 0.5460218 )+ 
    (ef_30_50 * 0.4191643) + 
    (ef_l30 * 1.094443) + 
    (bin_emergent  * 0.7127953) + 
    -4.789594)) %>% 
    mutate(euro_predmort = exp(z) / (1+exp(z)))
}

range_left = function(x, left, right){
  ifelse(x > left & x <= right, TRUE, FALSE)
}

# function to get time in range variable 
# inputs are df of hemodynamics, vector of thresholds, and hemodynamic variable (map or cvp)
get_ranges = function(hemo_data, thresholds, hemo_variable) {
  column_name = paste0("val_", hemo_variable) # get column of hemodynamic vars from hemo dataset
  
  hemo_data %>%
    # filter(cat_anes == "intra") %>%
    group_by(id) %>%
    select(id, time, all_of(column_name)) %>%
    mutate(across(
      all_of(column_name),
      ~ cut(.x,
            breaks = thresholds,
            include.lowest = TRUE)
    )) %>%
    count(across(all_of(column_name)), .drop = FALSE) %>%
    pivot_wider(
      values_from = n,
      id_cols = id,
      names_from = all_of(column_name)
    ) %>%
    ungroup() %>%
    rename(missing = 'NA') %>%
    rename_with(~ str_c(hemo_variable, "_", .),-id)
}


get_ranges_cpb = function(hemo_data, thresholds, hemo_variable) {
  column_name = paste0("val_", hemo_variable) # get column of hemodynamic vars from hemo dataset
  
  hemo_data %>%
    # filter(cat_anes == "intra") %>%
    group_by(id, cat_cpb) %>%
    select(id, time, cat_cpb, all_of(column_name)) %>%
    mutate(across(
      all_of(column_name),
      ~ cut(.x,
            breaks = thresholds,
            include.lowest = TRUE)
    )) %>%
    count(across(all_of(column_name)), .drop = FALSE) %>%
    pivot_wider(
      values_from = n,
      id_cols = c(id, cat_cpb),
      names_from = all_of(column_name)
    ) %>%
    ungroup() %>%
    rename(missing = 'NA') %>%
    rename_with(~ str_c(hemo_variable, "_", .),-c(id, cat_cpb))
}


get_ranges2 = function(hemo_data,
                          thresh1, var1,
                          thresh2, var2) {
  
  hemo_data %>%
    mutate(
      across({{var1}}, ~cut(.x, breaks = thresh1, include.lowest = TRUE), .names = "v1"),
      across({{var2}}, ~cut(.x, breaks = thresh2, include.lowest = TRUE), .names = "v2")) %>% 
    group_by(id) %>% 
    count(v1, v2, .drop = FALSE) %>% 
    drop_na() %>%
    ungroup() %>%
    mutate(var = paste0(sub(".*val\\_", "", var1), "_", v1, "_", sub(".*val\\_", "", var2), "_", v2)) %>% 
    ungroup() %>%
    pivot_wider(names_from = var,
                values_from = n,
                id_cols = id) %>%
    ungroup()
}

# function to perform regressions with output of get_ranges and covariate dataframe 
# inputs: 
# hemo_preds: hemodynamic variables from get_ranges function 
# covariate_df: covariates data frame
# covars_control: which covariates are we controlling for? 
# outcome: name of outcome in covariates dataframe, either bin_aki48h or bin_akicomposite
# adjustment_type: one of "multi" or "uni", i.e. adjustment for covariates or univariate (no adjustment) 
# cma: TRUE or FALSE - do correlation and multiplicity adjustment? 
# B: number of iterations for CMA 
# exponentiate: TRUE or FALSE - get results on OR scale? 
# analysis_units: effect size per xx minutes? default is 5
# hemo type: one of map or cvp

# helper function to fit model (multivariable) 

run_regression = function(hemo_preds,
                          covariate_df,
                          outcome,
                          covars_control,
                          adjustment_type = "multi",
                          cma = FALSE,
                          B = 1000,
                          exponentiate = TRUE,
                          analysis_units = 5,
                          hemo_type) {
  newnames = c(AKI = outcome) # to use select function
  hemo_vars = colnames(hemo_preds %>% select(-id)) # get hemodynamic variables
  # divide by analysis units for regression
  hemo_preds <- hemo_preds %>%
    mutate(across(
      starts_with(hemo_type) |
        starts_with("missing"),
      ~ .x / analysis_units
    ))
  fit_model = function(var, df) {
    data_temp = df %>%
      select(id, all_of(outcome), all_of(covars_control)) %>%
      rename(all_of(newnames)) %>%
      mutate(AKI = factor(AKI)) %>%
      left_join(hemo_preds %>% select(id, all_of(var)),
                by = join_by(id))
    
    logistic_reg() %>%
      set_engine("glm") %>%
      fit(AKI ~ .,
          data = data_temp %>% select(-id)) %>%
      broom::tidy() %>%
      filter(grepl(hemo_type, term))
  }
  
  # helper function to do univariate regression (no covariates)
  fit_model_uni = function(var, df) {
    data_temp = df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      mutate(AKI = factor(AKI)) %>%
      left_join(hemo_preds %>% select(id, all_of(var)),
                by = join_by(id))
    
    logistic_reg() %>%
      set_engine("glm") %>%
      fit(AKI ~ .,
          data = data_temp %>% select(-id)) %>%
      broom::tidy() %>%
      filter(grepl(hemo_type, term))
  }
  
  # perform "one at a time" regressions, store results
  if (adjustment_type == "uni") {
    summary = map_df(.x = hemo_vars,
                     .f = fit_model_uni,
                     df = covariate_df)
  }
  else{
    summary = map_df(.x = hemo_vars,
                     .f = fit_model,
                     df = covariate_df)
  }
  # do bonferroni correction
  n <- nrow(summary)
  bonferroni_p = 0.05 / n
  crit_val_bf <- qnorm(1 - (bonferroni_p / 2))
  
  
  # do fdr 
  fdr_pval = p.adjust(summary$p.value, method = "fdr")
  
  # and confidence intervals to summary
  summary <-
    summary %>%
    mutate(
      lb = estimate - (1.96 * std.error),
      ub = estimate + (1.96 * std.error),
      lb_bonf = estimate - (crit_val_bf * std.error),
      ub_bonf = estimate + (crit_val_bf * std.error),
      fdr_p = fdr_pval
    )
  # if cma is true, perform regression many times w replacement and get CMA critical val
  if (cma == TRUE) {
    aki_ids = covariate_df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      filter(AKI == 1) %>%
      pull(id)
    no_aki_ids = covariate_df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      filter(AKI == 0) %>%
      pull(id)
    
    sample_w_replacement <-
      function(vars) {
        samp_aki = sample(aki_ids, length(aki_ids), replace = TRUE)
        samp_noaki = sample(no_aki_ids, length(no_aki_ids), replace = TRUE)
        samp = c(samp_noaki, samp_aki)
        map_df(.x = vars,
               .f = fit_model,
               df = covariate_df %>% filter(id %in% samp)) %>%
          pull(estimate)
      }
    boostrap_coefs <-
      replicate(B,
                sample_w_replacement(hemo_vars)) %>%
      t()
    corr_tmp <- cor(boostrap_coefs)
    cma_val <- mvtnorm::qmvnorm(p = 0.95, corr = corr_tmp, tail = "both.tails")$quantile
    summary <-
      summary %>%
      mutate(
        lb_cma = estimate - (cma_val * std.error),
        ub_cma = estimate + (cma_val * std.error)
      )
    
  }
  if (exponentiate == TRUE) {
    summary <-
      summary %>%
      mutate(across(c(
        estimate, std.error, starts_with("lb"),
        starts_with("ub")
      ), exp))
  }
  summary
}

# function to get bivariate ranges (for map and cvp together) 
get_ranges_biv = function(hemo_data,
                          map_thresholds,
                          cvp_thresholds) {
  hemo_data %>%
    # filter(cat_anes == "intra") %>%
    group_by(id) %>%
    mutate(
      cut_map = cut(val_map, breaks = map_thresholds, include.lowest = TRUE),
      cut_cvp = cut(val_cvp, breaks = cvp_thresholds, include.lowest = TRUE)
    ) %>%
    count(cut_map, cut_cvp, .drop = F) %>%
    drop_na() %>%
    ungroup() %>%
    mutate(map_cvp = paste0("map", cut_map, "cvp", cut_cvp)) %>%
    ungroup() %>%
    pivot_wider(names_from = map_cvp,
                values_from = n,
                id_cols = id) %>%
    ungroup()
}

# function to run bivariate regression 
# same inputs as run_regression function except "hemo_type" not needed 

run_regression_biv = function(hemo_preds,
                              covariate_df,
                              outcome,
                              covars_control,
                              adjustment_type = "multi",
                              cma = FALSE,
                              B = 1000,
                              exponentiate = TRUE,
                              analysis_units = 5) {
  newnames = c(AKI = outcome) # for ease of select()
  hemo_vars = colnames(hemo_preds %>% select(-id)) # get hemodynamic variables
  # divide by analysis units for regression
  hemo_preds <- hemo_preds %>%
    mutate(across(!starts_with("id"), ~ .x / analysis_units))
  
  # helper function to fit regression with one variable at a time, controlling for covariates
  # extract just estimate for hemodynamic variable
  
  fit_model_bv = function(var, df) {
    data_temp = df %>%
      select(id, all_of(outcome), all_of(covars_control)) %>%
      rename(all_of(newnames)) %>%
      mutate(AKI = factor(AKI)) %>%
      left_join(hemo_preds %>% select(id, all_of(var)),
                by = join_by(id))
    
    logistic_reg() %>%
      set_engine("glm") %>%
      fit(AKI ~ .,
          data = data_temp %>% select(-id)) %>%
      broom::tidy() %>%
      filter(grepl("map", term))
  }
  # helper function to do univariate regression (no covariates)
  fit_model_uni_bv = function(var, df) {
    data_temp = df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      mutate(AKI = factor(AKI)) %>%
      left_join(hemo_preds %>% select(id, all_of(var)),
                by = join_by(id))
    
    logistic_reg() %>%
      set_engine("glm") %>%
      fit(AKI ~ .,
          data = data_temp %>% select(-id)) %>%
      broom::tidy() %>%
      filter(grepl("map", term))
  }
  # perform "one at a time" regressions, store results
  if (adjustment_type == "uni") {
    summary = map_df(.x = hemo_vars,
                     .f = fit_model_uni_bv,
                     df = covariate_df)
  }
  else{
    summary = map_df(.x = hemo_vars,
                     .f = fit_model_bv,
                     df = covariate_df)
  }
  # do bonferroni correction
  n <- nrow(summary)
  bonferroni_p = 0.05 / n
  crit_val <- qnorm(1 - (bonferroni_p/2))
  
  fdr_pval = p.adjust(summary$p.value, method = "fdr")
  
  # get confidence intervals 
  summary <-
    summary %>%
    mutate(
      lb = estimate - (1.96 * std.error),
      ub = estimate + (1.96 * std.error),
      lb_bonf = estimate - (crit_val * std.error),
      ub_bonf = estimate + (crit_val * std.error),
      fdr_p = fdr_pval
    )
  # if cma is true, perform regression many times w replacement and get CMA critical val
  if (cma == TRUE) {
    aki_ids = covariate_df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      filter(AKI == 1) %>%
      pull(id)
    no_aki_ids = covariate_df %>%
      select(id, all_of(outcome)) %>%
      rename(all_of(newnames)) %>%
      filter(AKI == 0) %>%
      pull(id)
    
    sample_w_replacement <-
      function(vars) {
        samp_aki = sample(aki_ids, length(aki_ids), replace = TRUE)
        samp_noaki = sample(no_aki_ids, length(no_aki_ids), replace = TRUE)
        samp = c(samp_noaki, samp_aki)
        map_df(.x = vars,
               .f = fit_model_bv,
               df = covariate_df %>% filter(id %in% samp)) %>%
          pull(estimate)
      }
    boostrap_coefs <-
      replicate(B,
                sample_w_replacement(hemo_vars)) %>%
      t()
    corr_tmp <- cor(boostrap_coefs)
    cma_val <- mvtnorm::qmvnorm(p = 0.95, corr = corr_tmp, tail = "both.tails")$quantile
    summary <-
      summary %>%
      mutate(
        lb_cma = estimate - (cma_val * std.error),
        ub_cma = estimate + (cma_val * std.error)
      )
    
  }
  if (exponentiate == TRUE) {
    summary <-
      summary %>%
      mutate(across(c(
        estimate, std.error, starts_with("lb"),
        starts_with("ub")
      ), exp))
  }
  summary
  
}


run_shells_univariate = function(var, df){
  tempdf = df %>% select(all_of(var), bin_aki48h)
  logistic_reg() %>%
    set_engine("glm") %>%
    fit(bin_aki48h ~ .,
        data = tempdf) %>% 
    tidy(exponentiate = TRUE) %>% 
    filter(grepl("group", term)) 
}
run_shells_multi = function(var, df){
  tempdf = df %>% select(id, all_of(var)) %>% 
    left_join(covar_data %>% select(id, bin_aki48h, all_of(primary_and_mediators))) %>% 
    select(-id) %>% 
    mutate(across(bin_aki48h, as.factor))
  logistic_reg() %>%
    set_engine("glm") %>%
    fit(bin_aki48h ~ .,
        data = tempdf) %>% 
    tidy(exponentiate = TRUE) %>% 
    filter(grepl("group", term)) 
  
}


