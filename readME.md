This repository contains the code for joint mapping of MAP and CI in the WAKE cohort of CABG surgeries and accompanies xx manuscript. 

Files are as follows: 

## code 

### Analysis pipeline 

+ `01_process_timeseries.R`
  + purpose: concatenate individual files into one and add time stamps
  + input: individual files for each subject w/ time series data from `data/raw/individual`
  + output: `data/processed/hemo_timeseries.rds`: one file w/ all time series data for all subjects
+ `02_interpolate.R`: 
  + purpose: filter based on rules and interpolate based on interpolation criteria
  + input: `data/processed/hemo_timeseries.rds`
  + output: `data/processed/hemo_timeseries_interp.rds`
+ `03_missingness.R`
  + purpose: figure out individuals to exclude for covariates or missingness criteria
  + input: `data/processed/hemo_timeseries_interp.rds`
  + output: `data/processed/exclusion_summary.rds`, `data/processed/hemo_analytic.rds`
+ `04_create_analysis_data.R`
  + purpose: calculate time in zones 
  + input: `data/processed/hemo_analytic.rds`
  + output: files in `data/analytic/` folder 
+ `05_process_covars.R`: 
  + purpose: process covariates
  + input: `data/raw/WAKE.flatfile.IDsfixed.11.27.24.xlsx`
  + output: `data/processed/covars_proc.rds`

### Utilities 

+ `utilities.R`: helpful functions for calculating time in ranges 
+ `filtering_settings.R`: parameters for filtering data 
+ `processing_fns.R`: functions for the data processing

### Manuscript

+ `manuscript_analyses_final.Rmd` and `manuscript_analyses_final.html`: all code for generating tables, figures, and results in manuscript 
