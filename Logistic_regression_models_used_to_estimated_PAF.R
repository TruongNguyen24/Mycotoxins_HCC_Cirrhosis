# Loading library

library(tidyverse)
library(sumvar)
library(splines)
library(survival)
library(survminer)
library(car)
library(rms)
library(haven)
library(ggdiagnostics)
library(graphPAF)
library(patchwork)
library(future.apply)
library(broom)

# Loading dataset

summary(af_merged) # This the tabular dataset of exposures: HBV, AFB1, HIV, HCV, smoking, tuberculosis, alcohol consumption, diabetes
af_merged$haz_alc <- ifelse(af_merged$who_audit_cat >= 2, 1, 0 )
af_merged$haz_alc <- factor(af_merged$haz_alc, levels = c(0,1))
af_merged$sex <- factor(af_merged$sex, levels = c(1,2))
af_merged$diabetes <- factor(af_merged$diabetes, levels = c(0,1))
af_merged$prevtb <- factor(af_merged$prevtb, levels = c(0,1))
af_merged$hiv_status <- factor(af_merged$hiv_status, levels = c(0,1))
af_merged$schisto_bin <- factor(af_merged$schisto_bin, levels = c(0,1))
af_merged$smoked_ever  <- factor(af_merged$smoked_ever , levels = c(0,1))
af_merged$hbsag_result   <- factor(af_merged$hbsag_result  , levels = c(0,1))
af_merged$hcvrna_result  <- factor(af_merged$hcvrna_result , levels = c(0,1))

# Setting environment

options(boot.parallel = "snow")
options(boot.ncpus = parallel::detectCores())

# Calculating weight for HBsAg
## 1) target prevalence from serosurvey
p_target <- 0.046 # prevalence of HBsAg positive in community serosurvey

## 2) observed prevalence among CONTROLS in your data
p_ctrl <- af_merged %>%
  filter(cc == "control") %>%
  dplyr::count(hbsag_result) %>%
  mutate(prop = n / sum(n)) %>%
  select(hbsag_result, prop) %>%   # <-- keep 2 columns
  tibble::deframe()                # named vector: c(Negative=..., Positive=...)

p_ctrl

## 3) compute post-strat weights for controls
w_pos <- p_target / p_ctrl["1"]
w_neg <- (1 - p_target) / p_ctrl["0"]

af_merged <- af_merged %>%
  mutate(
    hbsag_pos = hbsag_result %in% c("Positive","POS","pos", "1",TRUE),
    w = case_when(
      cc == "control" & hbsag_pos ~ w_pos,
      cc == "control" & !hbsag_pos ~ w_neg,
      TRUE ~ 1  # cases (or non-controls)
    )
  )

## Models for cirrhosis = cir_case

cc_cir <- af_merged[af_merged$pcat.x != 2, ] # Note: 1, Cirhosis; 2, HCC ; 3: control
cc_cir$cc <- factor(cc_cir$cc, levels = c("control", "Cirrhosis"))
cc_cir$cir_case <- cc_cir$cc

## Cirrhosis models
### Calculating ORs of cirrhosis (vs. control) for sex, age, HBV (HBsAg) 
lr_hbv_cirr <- glm(cir_case ~  final_result_HBsAg, 
              data = cc_cir, family = binomial)
tidy(lr_hbv_cirr, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = sprintf("%.6f", p.value))

lr_age_cirr <- glm(cir_case ~  age, 
                   data = cc_cir, family = binomial, weights = w)
tidy(lr_age_cirr, conf.int = T, exponentiate = T)

lr_sex_cirr <- glm(cir_case ~  sex, 
                   data = cc_cir, family = binomial, weights = w)
tidy(lr_sex_cirr, conf.int = T, exponentiate = T)

### Univariate analysis
#### AFB1

lr_as_af_w <- glm(cir_case ~  ns(age,2) + sex + Afb1_lysine_log, # AFB1 concentration is treated as continuous variable
                   data = cc_cir, family = binomial, weights = s_weights)
lr_as_af_w 
summary(lr_as_af_w) 

lr_as_af_w2 <- glm(cir_case ~  ns(age,2) + sex + af_bin, # AFB1 exposure is treated as binary variable
                   data = cc_cir, family = binomial, weights = w)
lr_as_af_w2 
summary(lr_as_af_w2)
tidy(lr_as_af_w2, conf.int = T, exponentiate = T)

system.time(
  af_cir <- PAF_calc_discrete(model=lr_as_af_w2, 
                               riskfactor= "af_bin", refval = 0,
                               data= cc_cir,
                               calculation_method = "B",
                               ci = TRUE, 
                               ci_type = "perc",
                               verbose= TRUE,
                               boot_rep = 5000)
)
print("1. Cirrhosis, aflatoxin")
print(af_cir) # Calculting PAF of AFB1 for Cirrhosis


#### HBV
lr_as_hbv <- glm(cir_case ~  ns(age,2) + sex +  hbsag_result,
                     data = cc_cir, family = binomial, weights = w)
tidy(lr_as_hbv, conf.int = T, exponentiate = T)

summary(lr_as_hbv)
system.time(
  hbv_cir <- PAF_calc_discrete(model=lr_as_hbv, 
                                 riskfactor= "hbsag_result", refval = 0,
                                 data= cc_cir,
                                 calculation_method = "B",
                                 ci = TRUE, 
                                 ci_type = "perc",
                                 verbose= TRUE,
                                 boot_rep = 5000)
)
print("2. Cirrhosis, HBV")
print(hbv_cir)

#### Smoking
lr_as_smoke_b <- glm(cir_case_bin ~  ns(age,2) + sex + smoke,
                     data = cc_cir, family = binomial, weights = w)

summary(lr_as_smoke_b)
tidy(lr_as_smoke_b, conf.int = T, exponentiate = T)

system.time(
  smoke_cir <- PAF_calc_discrete(model=lr_as_smoke_b, 
                                 riskfactor= "smoke", refval = 0,
                                 data= cc_cir,
                                 calculation_method = "B",
                                 ci = TRUE, 
                                 ci_type = "perc",
                                 verbose= TRUE,
                                 boot_rep = 5000)
)
print("2. Cirrhosis, smoking")
print(smoke_cir)

#### Diabetes
lr_as_dm <- glm(cir_case ~  ns(age,2) + sex + diabetes,
                data = cc_cir, family = binomial, weights = w)
summary(lr_as_dm)
tidy(lr_as_dm, conf.int = T, exponentiate = T)

system.time(
  dm_cir <- PAF_calc_discrete(model=lr_as_dm, 
                              riskfactor= "diabetes", refval = 0,
                              data= cc_cir,
                              calculation_method = "B",
                              ci = TRUE, 
                              ci_type = "perc",
                              verbose= TRUE,
                              boot_rep = 50)
)
print("3. Cirrhosis, diabetes")
print(dm_cir)

#### HIV
lr_as_hiv <- glm(cir_case ~  ns(age,2) + sex + hiv_status,
                 data = cc_cir, family = binomial, weights = w)
summary(lr_as_hiv)
tidy(lr_as_hiv, conf.int = T, exponentiate = T)

system.time(
  hiv_cir <- PAF_calc_discrete(model=lr_as_hiv, 
                               riskfactor= "hiv_status", refval = 0,
                               data= cc_cir,
                               calculation_method = "B",
                               ci = TRUE, 
                               ci_type = "perc",
                               verbose= TRUE,
                               boot_rep = 5000)
)
print("4. Cirrhosis, HIV")
print(hiv_cir)

#### Tuberculosis
lr_as_tb <- glm(cir_case ~  ns(age,2) + sex + prevtb,
                data = cc_cir, family = binomial, weights = w)
summary(lr_as_tb)
tidy(lr_as_tb, conf.int = T, exponentiate = T)

system.time(
  tb_cir <- PAF_calc_discrete(model=lr_as_tb, 
                                  riskfactor= "prevtb", refval = 0,
                                  data= cc_cir,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("6. Cirrhosis, TB")
print(tb_cir)

#### Hazardous alcohol consumption
lr_as_alc <- glm(cir_case ~  ns(age,2) + sex + haz_alc, # hazardous consumption
                data = cc_cir, family = binomial, weights = w)
summary(lr_as_alc)
tidy(lr_as_alc, conf.int = T, exponentiate = T)


lr_as_harm_alc <- glm(cir_case ~  ns(age,2) + sex + harm_alc, # harmful consumption
                data = cc_cir, family = binomial, weights = w)
summary(lr_as_harm_alc)

system.time(
  alc_cir <- PAF_calc_discrete(model=lr_as_alc, 
                                  riskfactor= "haz_alc", refval = 0,
                                  data= cc_cir,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("6. Cirrhosis, Hazardous alcohol")
print( alc_cir)

#### HCV
lr_as_hcv <- glm(cir_case ~  ns(age,2) + sex + hcvrna_result,
                data = cc_cir, family = binomial, weights = w)
summary(lr_as_hcv)
tidy(lr_as_hcv, conf.int = T, exponentiate = F)

system.time(
  hcv_cir <- PAF_calc_discrete(model=lr_as_hcv, 
                                  riskfactor= "hcvrna_result", refval = 0,
                                  data= cc_cir,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 10)
)
print("6. Cirrhosis,HCV")
print( hcv_cir)

### Multivariate analysis

lr_mv <- glm(cir_case ~  ns(age,2) + sex + hbsag_result,
             data = cc_cir %>% filter(!is.na(hiv_status)), family = binomial, weights = w)
tidy(lr_mv, conf.int = T, exponentiate = T)

lr_mv1 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv1, conf.int = T, exponentiate = T)

lr_mv2 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv2, conf.int = T, exponentiate = T)

lr_mv3 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + Afb1_lysine_log,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv3, conf.int = T, exponentiate = T)

lr_mv4 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + smoke,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv4, conf.int = T, exponentiate = T)

lr_mv5 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + smoke + Afb1_lysine_log,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv5, conf.int = T, exponentiate = T)

lr_mv6 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result + hiv_status + Afb1_lysine_log,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv6, conf.int = T, exponentiate = T)

lr_mv7 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv7, conf.int = T, exponentiate = T)

lr_mv7_lin <- glm(cir_case_bin ~  age + sex + hbsag_result * hiv_status + prevtb,
                  data = cc_cir, family = binomial, weights = w)
tidy(lr_mv7_lin, conf.int = T, exponentiate = T) 

lr_mv7_af <- glm(cir_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + af_bin,
              data = cc_cir, family = binomial, weights =  w)
tidy(lr_mv7_af, conf.int = T, exponentiate = T)

lr_mv8 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + diabetes,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv8, conf.int = T, exponentiate = T)

lr_mv9 <- glm(cir_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + smoked_ever,
              data = cc_cir, family = binomial, weights = w)
tidy(lr_mv9, conf.int = T, exponentiate = T)

#### Explore information criteira
AIC(lr_mv,
    lr_mv1,
    lr_mv2,
    lr_mv3,
    lr_mv4,
    lr_mv5,
    lr_mv6,
    lr_mv7_lin,
    lr_mv7,
    lr_mv7_af,
    lr_mv8,
    lr_mv9
    )
BIC(lr_mv,
    lr_mv1,
    lr_mv2,
    lr_mv3,
    lr_mv4,
    lr_mv5,
    lr_mv6,
    lr_mv7_lin,
    lr_mv7,
    lr_mv7_af,
    lr_mv8,
    lr_mv9)

#### Compare interaction terms with LRT
anova(lr_mv3, lr_mv7, test="LRT")
anova(lr_mv4, lr_mv9, test="LRT")

lr_mv7_af <- glm(cir_case ~  ns(age,2) + sex + hbsag_result* hiv_status + prevtb + af_bin,
              data = cc_cir, family = binomial, weights =  w)  # This is the final model for Cirrhosis

vif(lr_mv7_af)
summary(lr_mv7_af)
tidy(lr_mv7_af, conf.int = T, exponentiate = T) 


##### Estimate PAF for HBV in MV model
system.time(
  hbv_cir_tb <- PAF_calc_discrete(model=lr_mv7_af, 
                                  riskfactor= "hbsag_result", refval = 0,
                                  data= cc_cir,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  weight_vec = cc_cir$w,
                                  boot_rep = 5000)
)
print("1. MV Cirrhosis, HBV")
print(hbv_cir_tb)

##### Estimate PAF for AFB1 in MV model
system.time(
  PAF_af_cir <- PAF_calc_discrete(model=lr_mv7_af, 
                                  riskfactor= "af_bin", refval = 0,
                                  data= cc_cir,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  weight_vec = cc_cir$w,
                                  boot_rep = 5000)
)
print("1. MV Cirrhosis, HBV")
print( PAF_af_cir)

# Estimate joint/sequential PAF for HIV and TB (casually linked in DAG)
# Model list
 model_list <- list(
  glm(hiv_status ~ ns(age,2) + sex + hbsag_result +  af_bin, 
      data = cc_cir, family = binomial, weights = w),
  glm(prevtb  ~ ns(age,2) + sex + hbsag_result + hiv_status +  af_bin, 
      data = cc_cir, family = binomial, weights = w),
  glm(cir_case ~ ns(age,2) + sex + hbsag_result * hiv_status + prevtb + af_bin, 
      data = cc_cir, family = binomial, weights = w)
)

parent_list <- list(
  c("age","sex",  "hbsag_result", "af_bin"),
  c("age","sex","hbsag_result","hiv_status", "af_bin"),   
  c("age","sex","hbsag_result","hiv_status","prevtb", "af_bin")
)
node_vec <- c("hiv_status","prevtb","cir_case")

ap <- average_paf(
  data           = cc_cir,
  model_list     = model_list,
  parent_list    = parent_list,
  node_vec       = node_vec,
  riskfactor_vec = c("hiv_status","prevtb"), 
  ci             = TRUE,
  boot_rep       = 5000,
  ci_type        = "perc",
  ci_level       = 0.95,
  weight_vec     = cc_cir$w,
  verbose        = TRUE
)

print(ap)

plot(ap,max_PAF=0.5,min_PAF=-0.1,number_rows=3)

## Models for HCC = hcc_case

cc_hcc <- af_merged[af_merged$pcat.x != 1, ]
cc_hcc$cc <- factor(cc_hcc$cc, levels = c("control", "HCC"))
summary(cc_hcc)
cc_hcc$hcc_case <- cc_hcc$cc

## HCC models
### Calculating ORs of HCC (vs. control) for sex, age, HBV (HBsAg) 
lr_hbv_hcc <- glm(hcc_case ~  final_result_HBsAg, 
              data = cc_hcc, family = binomial)
tidy(lr_hbv_hcc, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = sprintf("%.6f", p.value))

lr_age_hcc <- glm(hcc_case ~  age, 
                   data = cc_hcc, family = binomial, weights = w)
tidy(lr_age_cirr, conf.int = T, exponentiate = T)

lr_sex_hcc <- glm(hcc_case ~  sex, 
                   data = cc_hcc, family = binomial, weights = w)
tidy(lr_sex_hcc, conf.int = T, exponentiate = T)

### Univariate analysis : similar to Cirrhosis univariate models

### Mutivariate analysis

lr_mv_hcc <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result,
             data = cc_hcc %>% filter(!is.na(hiv_status)), family = binomial, weights = w)
tidy(lr_mv_hcc, conf.int = T, exponentiate = T)

lr_mv_hcc_1 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_1, conf.int = T, exponentiate = T)

lr_mv_hcc_2 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_2, conf.int = T, exponentiate = T)

lr_mv_hcc_3 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + af_bin,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_3, conf.int = T, exponentiate = T)

lr_mv_hcc_4 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + smoke,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_4, conf.int = T, exponentiate = T)

lr_mv_hcc_5 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + prevtb + smoke + af_bin,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_5, conf.int = T, exponentiate = T)

lr_mv_hcc_6 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + af_bin,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_6, conf.int = T, exponentiate = T)

lr_mv_hcc_7 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_7, conf.int = T, exponentiate = T)

lr_mv7_hcc_lin <- glm(hcc_case ~  age + sex + hbsag_result * hiv_status + prevtb,
                  data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv7_hcc_lin, conf.int = T, exponentiate = T) 

lr_mv_hcc_7 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + af_bin,
              data = cc_hcc, family = binomial, weights =  w)
tidy(lr_mv_hcc_7, conf.int = T, exponentiate = T)

lr_mv_hcc_8 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + diabetes,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_8, conf.int = T, exponentiate = T)

lr_mv_hcc_9 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + smoked_ever,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_9, conf.int = T, exponentiate = T)

lr_mv_hcc_10 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result * hiv_status + prevtb + smoked_ever + af_bin,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_10, conf.int = T, exponentiate = T)

lr_mv_hcc_11 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + smoked_ever + haz_alc,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_11, conf.int = T, exponentiate = T)

lr_mv_hcc_12 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + smoked_ever + af_bin,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_12, conf.int = T, exponentiate = T)

lr_mv_hcc_13 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result + hiv_status + smoked_ever + af_bin + haz_alc,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_13, conf.int = T, exponentiate = T)

lr_mv_hcc_14 <- glm(hcc_case ~  ns(age,2) + sex + hbsag_result*hiv_status + smoke + haz_alc,
              data = cc_hcc, family = binomial, weights = w)
tidy(lr_mv_hcc_14, conf.int = T, exponentiate = T)

# Explore information criteira
AIC(lr_mv_hcc,
    lr_mv_hcc_1,
    lr_mv_hcc_2,
    lr_mv_hcc_3,
    lr_mv_hcc_4,
    lr_mv_hcc_5,
    lr_mv_hcc_6,
    lr_mv7_hcc_lin,
    lr_mv_hcc_7,
    lr_mv_hcc_8,
    lr_mv_hcc_9,
    lr_mv_hcc_10,
    lr_mv_hcc_11,
    lr_mv_hcc_12,
     lr_mv_hcc_13,
    lr_mv_hcc_14
    )
BIC(lr_mv_hcc,
    lr_mv_hcc_1,
    lr_mv_hcc_2,
    lr_mv_hcc_3,
    lr_mv_hcc_4,
    lr_mv_hcc_5,
    lr_mv_hcc_6,
    lr_mv7_hcc_lin,
    lr_mv_hcc_7,
    lr_mv_hcc_8,
    lr_mv_hcc_9,
    lr_mv_hcc_10,
    lr_mv_hcc_11,
    lr_mv_hcc_12,
    lr_mv_hcc_13,
    lr_mv_hcc_14)

# Compare interaction terms with LRT
anova(lr_mv_hcc_6, lr_mv_hcc_11, test="LRT")
anova(lr_mv_hcc_6, lr_mv_hcc_11, test="LRT")

# Final model selection
cc_hcc$hcc_case_bin <- factor(cc_hcc$hcc_case_bin, levels = c(0,1))

lr_mv_hcc_6 <- glm(hcc_case_bin ~  ns(age,2) + sex + hbsag_result + hiv_status + af_bin,
              data = cc_hcc, family = binomial, weights = s_weights)
tidy(lr_mv_hcc_6, conf.int = T, exponentiate = T)

summary(lr_mv7)
cc_cir$hbsag_result <- factor(cc_cir$hbsag_result, levels = c(0,1))
summary(cc_cir)

lr_mv_hcc_12 <- glm(hcc_case_bin ~  age + sex + hbsag_result + hiv_status + smoke + af_bin,
              data = cc_hcc, family = binomial, weights = w) # This is final model for HCC
tidy(lr_mv_hcc_12, conf.int = T, exponentiate = T)
vif(lr_mv_hcc_12)

system.time(
  hbv_hcc <- PAF_calc_discrete(model=lr_mv_hcc_12, 
                                  riskfactor= "hbsag_result", refval = 0,
                                  data= cc_hcc,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("1. MV HCC, HBV")
print(hbv_hcc)

system.time(
  hiv_hcc <- PAF_calc_discrete(model=lr_mv_hcc_12, 
                                  riskfactor= "hiv_status", refval = 0,
                                  data= cc_hcc,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("1. MV HCC, HIV")
print(hiv_hcc)

system.time(
  smoke_hcc <- PAF_calc_discrete(model=lr_mv_hcc_12, 
                                  riskfactor= "smoke", refval = 0,
                                  data= cc_hcc,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("1. MV HCC, smoke")
print(smoke_hcc)

system.time(
  af_hcc <- PAF_calc_discrete(model=lr_mv_hcc_12, 
                                  riskfactor= "af_bin", refval = 0,
                                  data= cc_hcc,
                                  calculation_method = "B",
                                  ci = TRUE, 
                                  ci_type = "perc",
                                  verbose= TRUE,
                                  boot_rep = 5000)
)
print("1. MV HCC, AFB1")
print(af_hcc)



