################################################################################
# Cost-effectiveness of leadless versus transvenous cardiac pacing
# in haemodialysis patients
#
# Public code version shared for transparency of the analytical workflow.
# Original SNDS/REIN individual-level data are confidential and cannot be shared.
#
# This script assumes that the secured analytic datasets have already been
# created in the authorized environment.
################################################################################

# ------------------------------------------------------------------------------
# REQUIRED ANALYTIC OBJECTS
# ------------------------------------------------------------------------------
# The code below assumes the following objects are already available:
#
# df_ps                  : analytic dataset used for propensity score estimation
# Table_finale_4         : main analytic dataset
# Table_finale_4_netben  : main analytic dataset with updated propensity score
# Table_finale_4_sdnt    : dataset with costs excluding dialysis/transportation
# Table_finale_4_hosp    : dataset with hospital costs only
# Table_cout_totaux_mensuel_323 : monthly cost categories dataset
# Cost_matrix            : total cost matrix
# Cost_matrix_sdnt       : cost matrix excluding dialysis/transportation
# Cost_matrix_hosp       : hospital-only cost matrix
# X_covars               : covariate matrix used in doubly robust NetBenReg models
# Part.times             : partition time vector for NetBenReg
# lambda_year            : vector of willingness-to-pay thresholds
# form_ps                : propensity score model formula
#
# Optional objects used for some figures:
# Table_effectiveness_2  : dataset including date_index
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------

library(tidyverse)
library(survival)
library(WeightIt)
library(cobalt)
library(NetBenReg)
library(survey)
library(ggplot2)
library(scales)

################################################################################
# UTILITY FUNCTIONS
################################################################################

get_INB_from_NetBenReg <- function(fit, lambda_target, trt_name = "group") {
  
  idx_lambda <- which(sapply(fit, function(x) {
    !is.null(x$Reg.type) &&
      x$Reg.type == "NBR" &&
      !is.null(x$lambda) &&
      !is.na(x$lambda) &&
      x$lambda == lambda_target
  }))
  
  if (length(idx_lambda) != 1) {
    stop(paste0("lambda = ", lambda_target, " not found uniquely"))
  }
  
  ct_lambda <- fit[[idx_lambda]]$coef.table
  
  est <- ct_lambda[trt_name, "Estimate"]
  se  <- ct_lambda[trt_name, "Std.err"]
  ci  <- est + c(-1, 1) * 1.96 * se
  
  p_ce <- if (is.na(se) || se <= 0) NA_real_ else pnorm(est / se)
  
  idx_eff <- which(sapply(fit, function(x) {
    !is.null(x$Reg.type) && x$Reg.type == "Effect"
  }))
  
  if (length(idx_eff) != 1) stop("Effect block not found uniquely")
  
  ct_eff <- fit[[idx_eff]]$coef.table
  
  delta_eff_est <- ct_eff[trt_name, "Estimate"]
  delta_eff_se  <- ct_eff[trt_name, "Std.err"]
  delta_eff_ci  <- delta_eff_est + c(-1, 1) * 1.96 * delta_eff_se
  
  if (lambda_target == 0) {
    delta_cost_est <- -est
    delta_cost_ci  <- -rev(ci)
    
    ICER_est <- if (!is.na(delta_eff_est) && abs(delta_eff_est) > 1e-6) {
      delta_cost_est / delta_eff_est
    } else {
      NA_real_
    }
  } else {
    delta_cost_est <- NA_real_
    delta_cost_ci  <- c(NA_real_, NA_real_)
    ICER_est       <- NA_real_
  }
  
  list(
    lambda             = lambda_target,
    INB_est            = est,
    INB_se             = se,
    INB_ci_low         = ci[1],
    INB_ci_high        = ci[2],
    p_cost_effective   = p_ce,
    delta_cost_est     = delta_cost_est,
    delta_cost_ci_low  = delta_cost_ci[1],
    delta_cost_ci_high = delta_cost_ci[2],
    delta_eff_est      = delta_eff_est,
    delta_eff_se       = delta_eff_se,
    delta_eff_ci_low   = delta_eff_ci[1],
    delta_eff_ci_high  = delta_eff_ci[2],
    ICER_est           = ICER_est
  )
}

extract_INB_multiple_lambda <- function(fit, lambdas, trt_name = "group") {
  do.call(rbind, lapply(lambdas, function(lam) {
    r <- get_INB_from_NetBenReg(fit, lambda_target = lam, trt_name = trt_name)
    data.frame(
      lambda = lam,
      INB_est = r$INB_est,
      INB_ci_low = r$INB_ci_low,
      INB_ci_high = r$INB_ci_high,
      p_cost_effective = r$p_cost_effective,
      delta_eff_est = r$delta_eff_est,
      row.names = NULL
    )
  }))
}

################################################################################
# A) PROPENSITY SCORE WEIGHTING
################################################################################

# ------------------------------------------------------------------------------
# PROPENSITY SCORE MODEL SPECIFICATION
# ------------------------------------------------------------------------------

# Covariates included in the propensity score model
ps_vars <- c(
  "age_inclusion",
  "SEX",
  "BMI",
  "duree_dial_cat",
  "MARCHn",
  "Cisch",
  "ACFA",
  "BAV",
  "HTAP",
  "TDRV",
  "TDRSV",
  "LIP",
  "TRICUSP",
  "BLOC",
  "DYSSIN",
  "nombre_PM_cat",
  "Abord_a_inclusion_etude"
)

# Categorical covariates among propensity score variables
ps_factor_vars <- c(
  "SEX",
  "duree_dial_cat",
  "MARCHn",
  "Cisch",
  "ACFA",
  "BAV",
  "HTAP",
  "TDRV",
  "TDRSV",
  "LIP",
  "TRICUSP",
  "BLOC",
  "DYSSIN",
  "nombre_PM_cat",
  "Abord_a_inclusion_etude"
)

# Additional descriptive variables displayed in balance diagnostics
extra_balance_vars <- c(
  "nephgp", "DIABn", "SASn", "ICn", "RYTHMn",
  "ANEVn", "AMIn", "AVCAITn", "KCn"
)

# Keep only variables available in the secured analytic dataset
ps_vars_in_df <- intersect(ps_vars, names(df))
ps_missing_vars <- setdiff(ps_vars, names(df))
if (length(ps_missing_vars) > 0) {
  warning("Variables removed from propensity score model: ",
          paste(ps_missing_vars, collapse = ", "))
}

ps_factor_vars_in_df <- intersect(ps_factor_vars, names(df))
extra_balance_vars_in_df <- intersect(extra_balance_vars, names(df))

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

# Convert categorical variables to factors
df <- df %>%
  mutate(across(all_of(ps_factor_vars_in_df), ~ factor(.))) %>%
  mutate(across(all_of(extra_balance_vars_in_df), ~ factor(.)))

# Separate continuous PS variables
ps_numeric_vars_in_df <- setdiff(ps_vars_in_df, ps_factor_vars_in_df)

# Working copy for propensity score estimation
df_ps <- df

# Median imputation for continuous covariates
for (v in ps_numeric_vars_in_df) {
  med_v <- median(df_ps[[v]], na.rm = TRUE)
  df_ps[[v]][is.na(df_ps[[v]])] <- med_v
}

# Explicit "Missing" level for categorical variables
add_missing_level <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "Missing"
  factor(x)
}

df_ps <- df_ps %>%
  mutate(across(all_of(ps_factor_vars_in_df), add_missing_level)) %>%
  mutate(across(all_of(extra_balance_vars_in_df), add_missing_level))

# Treatment indicator as factor
df_ps <- df_ps %>%
  mutate(
    treated = factor(
      treated,
      levels = c(0, 1),
      labels = c("Transvenous pacemaker", "Micra")
    )
  )

# ------------------------------------------------------------------------------
# PROPENSITY SCORE ESTIMATION
# ------------------------------------------------------------------------------

form_ps <- as.formula(
  paste("treated ~", paste(ps_vars_in_df, collapse = " + "))
)

w.out <- weightit(
  formula  = form_ps,
  data     = df_ps,
  method   = "ps",
  estimand = "ATT"
)

summary(w.out)

# Active propensity score used in the analyses
df_ps$propensity_score <- w.out$ps

# Quick checks
summary(df_ps$propensity_score)
hist(df_ps$propensity_score)

# ------------------------------------------------------------------------------
# COVARIATE BALANCE
# ------------------------------------------------------------------------------

addl_form <- if (length(extra_balance_vars_in_df) > 0) {
  as.formula(paste("~", paste(extra_balance_vars_in_df, collapse = " + ")))
} else {
  NULL
}

bal <- bal.tab(
  x           = w.out,
  data        = df_ps,
  un          = TRUE,
  disp.means  = TRUE,
  m.threshold = 0.15,
  stats       = "mean.diffs",
  binary      = "std",
  abs         = TRUE,
  addl        = addl_form
)

bal

var_names <- c(
  age_inclusion = "Age at inclusion",
  SEX_M = "Male sex",
  BMI = "Body mass index",

  `duree_dial_cat_[   2, 290)` = "Dialysis duration: 2–290 days",
  `duree_dial_cat_[ 290, 775)` = "Dialysis duration: 290–775 days",
  `duree_dial_cat_[ 775,1846)` = "Dialysis duration: 775–1846 days",
  `duree_dial_cat_[1846,4811]` = "Dialysis duration: ≥1846 days",

  `MARCHn_Difficulté march` = "Walking difficulty",
  MARCHn_Marche_autonome = "Autonomous walking",
  MARCHn_Missing = "Walking ability: Missing",

  Cisch = "Ischemic heart disease",
  ACFA = "Atrial fibrillation",
  BAV = "Atrioventricular block",
  HTAP = "Pulmonary hypertension",
  TDRV = "Right ventricular dysfunction",
  TDRSV = "Left ventricular dysfunction",
  LIP = "Dyslipidemia",
  TRICUSP = "Tricuspid valve disease",
  BLOC = "Intraventricular conduction disorder",
  DYSSIN = "Cardiac dyssynchrony",

  AMIn = "History of myocardial infarction",
  DIABn = "Diabetes",
  SASn = "Obstructive sleep apnea",
  AVCAITn = "History of stroke or transient ischemic attack",
  ANEVn = "Aneurysm",
  ICn = "Heart failure",
  KCn = "History of cancer",

  RYTHMn_0 = "Rhythm disorder: No",
  RYTHMn_1 = "Rhythm disorder: Yes",
  RYTHMn_Missing = "Rhythm disorder: Missing",

  `nombre_PM_cat_[ 2, 28)` = "Center pacemaker implantation volume: 2–28",
  `nombre_PM_cat_[28, 48)` = "Center pacemaker implantation volume: 28–48",
  `nombre_PM_cat_[48, 78)` = "Center pacemaker implantation volume: 48–78",
  `nombre_PM_cat_[78,362]` = "Center pacemaker implantation volume: ≥78",

  Abord_a_inclusion_etude_Arteriovenous_fistula = "Arteriovenous fistula",
  Abord_a_inclusion_etude_Catheter = "Catheter",
  Abord_a_inclusion_etude_Other = "Other access",

  `nephgp_Glomérulonéphritis` = "Primary kidney disease: Glomerulonephritis",
  `nephgp_Hypertension` = "Primary kidney disease: Hypertension",
  `nephgp_Diabètes` = "Primary kidney disease: Diabetes",
  `nephgp_Autre` = "Primary kidney disease: Other",
  `nephgp_Inconnu` = "Primary kidney disease: Unknown",
  `nephgp_Polykystose` = "Primary kidney disease: Polycystic kidney disease"
)

vars_lp <- rownames(bal$Balance)
vars_lp <- vars_lp[!grepl("Missing", vars_lp, fixed = TRUE)]
vars_lp <- vars_lp[vars_lp != "MARCHn:<NA>"]

love.plot(
  bal,
  var.names    = var_names,
  var.include  = vars_lp,
  abs          = TRUE,
  binary       = "std",
  thresholds   = c(m = 0.15),
  colors       = c("red", "blue"),
  stars        = "raw",
  var.order    = "unadjusted",
  sample.names = c("Unweighted", "Weighted"),
  stat         = "mean.diffs",
  title        = "Standardized mean differences before and after weighting",
  xlab         = "|Standardized Mean Difference|"
)

# Weighted survey design object
des_iptw <- svydesign(
  ids = ~1,
  weights = ~w.out$weights,
  data = df_ps
)

# ------------------------------------------------------------------------------
# EXPORT PROPENSITY SCORE TO THE MAIN ANALYTIC DATASET
# ------------------------------------------------------------------------------

ps_export <- df_ps %>%
  dplyr::select(NUM_ENQ, propensity_score) %>%
  dplyr::distinct(NUM_ENQ, .keep_all = TRUE)

Table_finale_4_netben <- Table_finale_4 %>%
  dplyr::select(-dplyr::any_of("propensity_score")) %>%
  dplyr::left_join(ps_export, by = "NUM_ENQ") %>%
  mutate(
    treated_num = ifelse(treated %in% c(1, "Micra"), 1, 0)
  )

# ------------------------------------------------------------------------------
# COVARIATE MATRIX FOR DOUBLY ROBUST NETBENREG MODELS
# ------------------------------------------------------------------------------

covars_df <- df_ps %>%
  mutate(
    treated_num = ifelse(treated == "Micra", 1, 0)
  )

X_covars <- model.matrix(form_ps, data = covars_df)[, -1, drop = FALSE]

dim(X_covars)

################################################################################
# B) MAIN MEDICO-ECONOMIC MODELS: OVERALL SURVIVAL
################################################################################

# Total costs
fit_PT_adj_DR_SP <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_SW_adj_DR_SP <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

# Costs excluding dialysis and transportation
fit_PT_adj_DR_SP_sdnt <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix_sdnt,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_SW_adj_DR_SP_sdnt <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix_sdnt,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

# Hospital costs only
fit_PT_adj_DR_SP_hosp <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix_hosp,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_SW_adj_DR_SP_hosp <- NetBenReg(
  Followup      = Table_finale_4_netben$Followup_years,
  delta         = Table_finale_4_netben$DECES,
  group         = Table_finale_4_netben$treated_num,
  Cost          = Cost_matrix_hosp,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Z             = X_covars,
  PS.Z          = X_covars,
  Doubly.Robust = TRUE,
  lambda        = lambda_year,
  L             = max(Part.times)
)

################################################################################
# C) DEVICE-RELATED INFECTION MODELS
################################################################################

fit_endo_SW_DR <- NetBenReg(
  Followup      = Table_finale_4$delay_endocarditis_years,
  delta         = Table_finale_4$endocarditis2025_ME2,
  group         = Table_finale_4$treated,
  Cost          = Cost_matrix,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_endo_SW_DR_sdnt <- NetBenReg(
  Followup      = Table_finale_4_sdnt$delay_endocarditis_years,
  delta         = Table_finale_4_sdnt$endocarditis2025_ME2,
  group         = Table_finale_4_sdnt$treated,
  Cost          = Cost_matrix_sdnt,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_endo_SW_DR_hosp <- NetBenReg(
  Followup      = Table_finale_4_hosp$delay_endocarditis_years,
  delta         = Table_finale_4_hosp$endocarditis2025_ME2,
  group         = Table_finale_4_hosp$treated,
  Cost          = Cost_matrix_hosp,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "SW",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_endo_PT_DR <- NetBenReg(
  Followup      = Table_finale_4$delay_endocarditis_years,
  delta         = Table_finale_4$endocarditis2025_ME2,
  group         = Table_finale_4$treated,
  Cost          = Cost_matrix,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_endo_PT_DR_sdnt <- NetBenReg(
  Followup      = Table_finale_4_sdnt$delay_endocarditis_years,
  delta         = Table_finale_4_sdnt$endocarditis2025_ME2,
  group         = Table_finale_4_sdnt$treated,
  Cost          = Cost_matrix_sdnt,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

fit_endo_PT_DR_hosp <- NetBenReg(
  Followup      = Table_finale_4_hosp$delay_endocarditis_years,
  delta         = Table_finale_4_hosp$endocarditis2025_ME2,
  group         = Table_finale_4_hosp$treated,
  Cost          = Cost_matrix_hosp,
  Eff           = NULL,
  Eff.only      = TRUE,
  Part.times    = Part.times,
  Method        = "PT",
  Doubly.Robust = TRUE,
  Z             = X_covars,
  PS.Z          = X_covars,
  lambda        = lambda_year,
  L             = max(Part.times)
)

################################################################################
# D) EXTRACTION OF MAIN RESULTS
################################################################################

# Overall survival, lambda = 0
res_PT_total_0  <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP, 0)
res_SW_total_0  <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP, 0)
res_PT_sdnt_0   <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP_sdnt, 0)
res_SW_sdnt_0   <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP_sdnt, 0)
res_PT_hosp_0   <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP_hosp, 0)
res_SW_hosp_0   <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP_hosp, 0)

# Overall survival, lambda = 50,000
res_PT_total_50k <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP, 50000)
res_SW_total_50k <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP, 50000)
res_PT_sdnt_50k  <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP_sdnt, 50000)
res_SW_sdnt_50k  <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP_sdnt, 50000)
res_PT_hosp_50k  <- get_INB_from_NetBenReg(fit_PT_adj_DR_SP_hosp, 50000)
res_SW_hosp_50k  <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP_hosp, 50000)

# Device-related infection, lambda = 0
res_endo_SW_tot_0  <- get_INB_from_NetBenReg(fit_endo_SW_DR, 0)
res_endo_SW_sdnt_0 <- get_INB_from_NetBenReg(fit_endo_SW_DR_sdnt, 0)
res_endo_SW_hosp_0 <- get_INB_from_NetBenReg(fit_endo_SW_DR_hosp, 0)

res_endo_PT_tot_0  <- get_INB_from_NetBenReg(fit_endo_PT_DR, 0)
res_endo_PT_sdnt_0 <- get_INB_from_NetBenReg(fit_endo_PT_DR_sdnt, 0)
res_endo_PT_hosp_0 <- get_INB_from_NetBenReg(fit_endo_PT_DR_hosp, 0)

# Device-related infection, lambda = 50,000
res_endo_SW_tot_50k  <- get_INB_from_NetBenReg(fit_endo_SW_DR, 50000)
res_endo_SW_sdnt_50k <- get_INB_from_NetBenReg(fit_endo_SW_DR_sdnt, 50000)
res_endo_SW_hosp_50k <- get_INB_from_NetBenReg(fit_endo_SW_DR_hosp, 50000)

res_endo_PT_tot_50k  <- get_INB_from_NetBenReg(fit_endo_PT_DR, 50000)
res_endo_PT_sdnt_50k <- get_INB_from_NetBenReg(fit_endo_PT_DR_sdnt, 50000)
res_endo_PT_hosp_50k <- get_INB_from_NetBenReg(fit_endo_PT_DR_hosp, 50000)

################################################################################
# E) WEIGHTED COST DESCRIPTION BY CATEGORY
################################################################################

dialysis_cols <- c(
  "HD CENTRE", "HD ENT", "HD CENTRE-UDM",
  "UDM", "UDM ENT",
  "DIALYSE TIERCE",
  "DP", "DPCA", "DPCA ENT", "AUTODIALYSE ASS",
  "DPA", "MATERIEL"
)

hospital_cols <- c("HOSPITALISATION")
procedure_cols <- c("ACTES TECHNIQUE")
biology_cols <- c("BIOLOGIE ANAPAT")
paramedical_cols <- c("PARAMEDICAUX")
transport_cols <- c("TRANSPORT")
pharmacy_cols <- c("PHARMACIE")
other_cols <- c(
  "TELESURVEILLANC", "PRESTATION FINA",
  "Z - AUTRE", "ZZZZ", "REMUNERATION ME"
)

cost_patient_cat <- Table_cout_totaux_mensuel_323 %>%
  group_by(NUM_ENQ) %>%
  summarise(
    cout_dialysis = sum(rowSums(across(all_of(dialysis_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_hospital = sum(rowSums(across(all_of(hospital_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_procedures = sum(rowSums(across(all_of(procedure_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_biology = sum(rowSums(across(all_of(biology_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_paramedical = sum(rowSums(across(all_of(paramedical_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_transport = sum(rowSums(across(all_of(transport_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_pharmacy = sum(rowSums(across(all_of(pharmacy_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_other = sum(rowSums(across(all_of(other_cols)), na.rm = TRUE), na.rm = TRUE),
    cout_total = sum(montant_mensuel_total, na.rm = TRUE),
    .groups = "drop"
  )

cost_patient_cat_2 <- Table_finale_4_netben %>%
  semi_join(cost_patient_cat, by = "NUM_ENQ") %>%
  left_join(cost_patient_cat, by = "NUM_ENQ")

p_treated <- mean(cost_patient_cat_2$treated_num, na.rm = TRUE)

cost_patient_cat_2 <- cost_patient_cat_2 %>%
  mutate(
    weight = case_when(
      treated_num == 1 ~  p_treated / propensity_score,
      treated_num == 0 ~ (1 - p_treated) / (1 - propensity_score)
    )
  )

costs_weighted_summary <- cost_patient_cat_2 %>%
  group_by(treated_num) %>%
  summarise(
    tot_w    = sum(weight, na.rm = TRUE),
    tot_w_fu = sum(weight * Followup_years, na.rm = TRUE),
    across(
      starts_with("cout_"),
      list(
        w_mean = ~ weighted.mean(.x, w = weight, na.rm = TRUE),
        rate   = ~ sum(weight * .x, na.rm = TRUE) / tot_w_fu
      ),
      .names = "{fn}_{.col}"
    ),
    mean_fu_w = tot_w_fu / tot_w,
    .groups = "drop"
  )

plot_costs <- costs_weighted_summary %>%
  mutate(
    treated = factor(
      treated_num,
      levels = c(0, 1),
      labels = c("Transvenous pacemaker", "Leadless pacemaker")
    )
  ) %>%
  dplyr::select(treated, starts_with("w_mean_cout_")) %>%
  pivot_longer(
    cols      = dplyr::starts_with("w_mean_cout_"),
    names_to  = "category",
    values_to = "mean_cost"
  ) %>%
  mutate(
    category = gsub("^w_mean_cout_", "", category),
    category = factor(
      category,
      levels = c("total", "hospital", "dialysis", "procedures",
                 "biology", "pharmacy", "paramedical",
                 "transport", "other"),
      labels = c("Total", "Hospitalizations", "Dialysis", "Procedures",
                 "Laboratory/Pathology", "Pharmacy", "Paramedical care",
                 "Transportation", "Other")
    )
  )

ggplot(plot_costs, aes(x = category, y = mean_cost, fill = treated)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c(
      "Transvenous pacemaker" = "#4A86C5",
      "Leadless pacemaker"    = "#1F4E79"
    )
  ) +
  labs(
    x = "Cost category",
    y = "Weighted mean cost per patient (€)",
    fill = "Device",
    title = "Weighted mean costs per patient by category and pacemaker type"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

################################################################################
# F) CUMULATIVE WEIGHTED COST CURVES (COSTS EXCLUDING DIALYSIS/TRANSPORTATION)
################################################################################

iptw_key <- df_ps %>%
  dplyr::select(NUM_ENQ, treated) %>%
  dplyr::distinct(NUM_ENQ, .keep_all = TRUE) %>%
  mutate(w = w.out$weights) %>%
  transmute(NUM_ENQ, treated, w)

iptw_key <- iptw_key %>%
  mutate(
    treated = factor(
      treated,
      levels = c("Transvenous pacemaker", "Micra")
    )
  )

dialysis_excluded_cols <- c(
  "HD CENTRE", "HD ENT", "HD CENTRE-UDM",
  "UDM", "UDM ENT",
  "DIALYSE TIERCE",
  "DP", "DPCA", "DPCA ENT", "AUTODIALYSE ASS",
  "DPA", "MATERIEL"
)

transport_excluded_cols <- c("TRANSPORT")

cat_cols <- setdiff(
  names(Table_cout_totaux_mensuel_323),
  c("NUM_ENQ", "YYYYMM")
)

keep_cols <- setdiff(cat_cols, c(dialysis_excluded_cols, transport_excluded_cols))

cost_long_sdnt <- Table_cout_totaux_mensuel_323 %>%
  semi_join(iptw_key, by = "NUM_ENQ") %>%
  mutate(
    cost_sdnt = rowSums(across(all_of(keep_cols)), na.rm = TRUE),
    date_mois = as.Date(paste0(YYYYMM, "01"), format = "%Y%m%d")
  ) %>%
  left_join(
    Table_finale_4 %>% dplyr::select(NUM_ENQ, date_index, Followup_years),
    by = "NUM_ENQ"
  ) %>%
  mutate(
    mois_suivi = floor(
      as.numeric(difftime(date_mois, date_index, units = "days")) / 30.44
    ) + 1,
    followup_months = floor(Followup_years * 12)
  ) %>%
  filter(
    mois_suivi >= 1,
    mois_suivi <= 48,
    mois_suivi <= followup_months
  ) %>%
  group_by(NUM_ENQ, mois_suivi) %>%
  summarise(
    cost_sdnt = sum(cost_sdnt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(iptw_key, by = "NUM_ENQ") %>%
  mutate(
    weighted_cost = cost_sdnt * w
  )

plot_df <- cost_long_sdnt %>%
  group_by(treated, mois_suivi) %>%
  summarise(
    sum_w  = sum(w, na.rm = TRUE),
    sum_wc = sum(weighted_cost, na.rm = TRUE),
    mean_month = sum_wc / sum_w,
    .groups = "drop"
  ) %>%
  arrange(treated, mois_suivi) %>%
  group_by(treated) %>%
  mutate(
    mean_cum = cumsum(mean_month)
  ) %>%
  ungroup()

ggplot(plot_df, aes(x = mois_suivi, y = mean_cum, group = treated, linetype = treated)) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "Months since implantation",
    y = "Cumulative ATT-weighted mean cost (€)",
    title = "Cumulative monthly costs among patients still under follow-up (excluding dialysis and transportation)",
    linetype = "Device"
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top"
  )

################################################################################
# G) CEAC FIGURES
################################################################################

res_50k_surv <- get_INB_from_NetBenReg(fit_SW_adj_DR_SP_sdnt, 50000)
p_ref_surv <- res_50k_surv$p_cost_effective

res_50k_endo <- get_INB_from_NetBenReg(fit_endo_SW_DR_sdnt, 50000)
p_ref_endo <- res_50k_endo$p_cost_effective

graphics.off()
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

plot(
  fit_SW_adj_DR_SP_sdnt,
  main = "A. Overall survival",
  xlab = "Willingness-to-pay threshold (€)",
  ylab = "Probability of being cost-effective"
)
abline(v = 50000, lty = 2, lwd = 2)
abline(h = p_ref_surv, lty = 2, lwd = 2)
points(50000, p_ref_surv, pch = 19)

plot(
  fit_endo_SW_DR_sdnt,
  main = "B. Device-related infection-free survival",
  xlab = "Willingness-to-pay threshold (€)",
  ylab = "Probability of being cost-effective"
)
abline(v = 50000, lty = 2, lwd = 2)
abline(h = p_ref_endo, lty = 2, lwd = 2)
points(50000, p_ref_endo, pch = 19)

################################################################################
# H) INCREMENTAL NET BENEFIT CURVE
################################################################################

inb_df <- lapply(fit_SW_adj_DR_SP_sdnt, function(x) {
  if (is.null(x$Reg.type) || x$Reg.type != "NBR") return(NULL)
  
  ct <- x$coef.table
  
  data.frame(
    lambda = x$lambda,
    est    = ct["group", "Estimate"],
    se     = ct["group", "Std.err"]
  )
}) %>%
  bind_rows() %>%
  mutate(
    lcl = est - 1.96 * se,
    ucl = est + 1.96 * se
  ) %>%
  arrange(lambda)

ggplot(inb_df, aes(x = lambda, y = est)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Willingness-to-pay threshold (€/life-year)",
    y = "Incremental net benefit – Leadless vs transvenous pacemaker",
    title = "Incremental net benefit curve: primary analysis (costs excluding dialysis and transportation)"
  ) +
  theme_minimal()
