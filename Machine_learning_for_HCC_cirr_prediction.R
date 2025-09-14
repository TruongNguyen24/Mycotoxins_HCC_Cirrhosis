# Loading library
library(mlr)
library(limma)
library(dplyr)
library(Matrix)        # for xgboost
library(ROSE)
library(DMwR)
library(archdata) # archdata: Example Datasets from Archaeological Research
library(doParallel) # doParallel: Foreach Parallel Adaptor for the 'parallel' Package
library(Boruta) # Boruta: Wrapper Algorithm for All Relevant Feature Selection
library(broom) # broom: Convert Statistical Objects into Tidy Tibbles
library(caret) # caret: Classification and Regression Training
library(Ckmeans.1d.dp) # Ckmeans.1d.dp : Optimal, Fast, and Reproducible Univariate Clustering
library(ComplexHeatmap) # ComplexHeatmap: Making Complex Heatmaps
library(corrr) # corrr: Correlations in R
library(dplyr) # dplyr: A Grammar of Data manipulation
library(ggraph) # ggraph: An Implementation of Grammar of Graphics for Graphs and Networks
library(haven) # haven: Import and Export 'SPSS', 'Stata' and 'SAS' Files
#library(igraph) # igraph: Network Analysis and Visualization
library(knitr) # knitr: A General-Purpose Package for Dynamic Report Generation in R
# library(nVennR) # nVennR: Create n-Dimensional, Quasi-Proportional Venn Diagrams
library(missForest) # missForest: Nonparametric Missing Value Imputation using Random Forest
library(openxlsx) # openxlsx: Read, Write and Edit xlsx Files
library(pROC) # pROC: Display and Analyze ROC Curves
library(purrr) # purrr: Functional Programming Tools
library(readxl) # readxl: Read Excel Files
library(tibble) # tibble: Simple Data Frames
library(tidyverse) # tidyverse: Easily Install and Load the 'Tidyverse'
#install.packages("xgboost")
library(xgboost) # xgboost: Extreme Gradient Boosting
library(randomForest)# randomForest: Breiman and Cutler's Random Forests for Classification and Regression
#library(kknn)
#install.packages("kernlab")
library(kernlab)
#install.packages("e1071")
library(e1071)
#install.packages("plotROC")

# Preparing data
se_50
meta_mat <- t(assay(se_50,"logcounts"))

meta_df <- as.data.frame(meta_mat)
head(meta_df)  
colnames(meta_df) <- paste0("meta_", colnames(meta_df))

clin_df <- as.data.frame(colData(se_50)) |>
          dplyr::select(sex, age, bmi_new, Afb1_lysine_log, CIT_log, OTA_log, liver_control,hbsag_result, smoked_ever, run_id) 
table(is.na(clin_df$smoked_ever)) 


dim(clin_df)

stopifnot(identical(rownames(clin_df), rownames(meta_df)))  # same order?
data <- cbind(clin_df, meta_df)
dim(data)
summary(data)

data$smoked_ever <- as.factor(data$smoked_ever)
data$smoked_ever[is.na(data$smoked_ever)] <- 0 
data$run_id <- as.factor(data$run_id)

# make target a factor and declare "disease" as positive class
data$disease <- ifelse(data$liver_control == "Liver diseases", 1, 0)
data$disease <- factor(data$disease,
                       levels = c(0, 1),
                       labels = c("control", "liver diseases"))


factor_cols <- c("sex", "hbsag_result", "smoked_ever")       # exclude target

data <- createDummyFeatures(
          data,
          target = "disease",
          cols   = factor_cols,
          method = "reference"   # <-- keeps K-1 columns
        )

data <- data %>% select(-liver_control)
head(data)


# ------------ 2.  variable groups ----------------------------
## edit as needed
clinical_vars <- c("sex.Female", "age", "bmi_new", "Afb1_lysine_log", "CIT_log", "OTA_log", "hbsag_result.1", "smoked_ever.1")

meta_prefix   <- "meta_"                               # all metabolite peaks start with this

#######################################################################
#  HCC vs Cirrhosis vs Control – limma feature-selection + ML (mlr v2) #
#######################################################################

## ------------ 0. USER SETTINGS --------------------------------------
se_50                                 # SummarizedExperiment
clinical_vars <- c("sex", "age", "bmi_new",
                   "Afb1_lysine_log", "CIT_log", "OTA_log",
                   "hbsag_result", "smoked_ever")
target_var  <- "diagnosis"            # ▲ three-level factor
batch_var   <- "run_id"
meta_prefix <- "meta_"
top_n       <- 40
set.seed(42)

## ------------ 1. LIBRARIES ------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(mlr)            # v2.x
  library(dplyr)
  library(Matrix)         # for xgboost
})

## ------------ 2. LOAD + FLATTEN DATA --------------------------------
meta <- t(assay(se_50, 2)) |> as.data.frame()
colnames(meta) <- paste0(meta_prefix, colnames(meta))

clin <- as.data.frame(colData(se_50)) 
clin$diagnosis <- ifelse(clin$pcat == 1, 1, ifelse(clin$pcat == 2, 2, 3))
table(clin$diagnosis)

clin$diagnosis <- factor(clin$diagnosis, levels = c(1,2,3), labels = c("HCC", "Cirrhosis", "Control"))

clin <- clin |> select(all_of(c(clinical_vars, batch_var, target_var)))

stopifnot(identical(rownames(clin), rownames(meta)))
data <- cbind(clin, meta)

## ---- 2a. keep batch vector, then DROP as predictor -----------------
batch <- factor(data[[batch_var]])
names(batch) <- rownames(data)
data  <- data %>% select(-all_of(batch_var))

## ---- 2b. logical → factor, one-hot (K–1) ---------------------------
data <- data %>% mutate(across(where(is.logical),
                               ~ factor(ifelse(.x, "yes", "no"),
                                        levels = c("no", "yes"))))

fac_cols <- names(data)[sapply(data, is.factor) & names(data) != target_var]
data <- createDummyFeatures(data, target = target_var,
                            cols = fac_cols, method = "reference")
data$smoked_ever.1[is.na(data$smoked_ever.1)] <- 0

## ---- 2c. finalise target ------------------------------------------
data[[target_var]] <- factor(data[[target_var]],
                             levels = c("Control", "Cirrhosis", "HCC"))   # ▲

## ------------ 3. WRAPPER (limma → impute → scale) -------------------
makeLimmaWrapper <- function(learner) {

  trainfun <- function(data, target, args) {

    ## 3.1 metabolite columns & limma design
    meta_cols <- grep(paste0("^", meta_prefix), colnames(data), value = TRUE)
    X         <- t(as.matrix(data[, meta_cols, drop = FALSE]))    # features × samples

    grp       <- factor(data[[target]],
                        levels = c("Control", "Cirrhosis", "HCC"))        # ▲
    design    <- model.matrix(~ grp + batch[rownames(data)])              # ▲

    fit <- eBayes(lmFit(X, design))
    ## rank by overall F-test across the 3 groups
    tt  <- topTable(fit, coef = 2:3, number = Inf,
                    adjust.method = "BH", sort.by = "F")           # ▲
    tt$score <- tt$F                                                # ▲
    top_meta <- rownames(tt)[seq_len(top_n)]

    keep <- intersect(c(clinical_vars, top_meta), colnames(data))
    df   <- data[, c(target_var, keep), drop = FALSE]

    ## impute centre/scale numeric predictors
    num <- sapply(df, is.numeric)
    ctr <- colMeans(df[, num, drop = FALSE])
    scl <- apply(df[, num, drop = FALSE], 2, sd)
    scl[scl == 0] <- 1
    df[, num] <- sweep(df[, num, drop = FALSE], 2, ctr, "-")
    df[, num] <- sweep(df[, num, drop = FALSE], 2, scl, "/")

    list(data = df,
         control = list(keep = keep, num = num,
                        ctr = ctr, scl = scl))
  }

  predictfun <- function(data, target, args, control) {

    keep <- intersect(control$keep, colnames(data))
    new  <- data[ ,  keep, drop = FALSE]

    if (any(control$num)) {
      nms <- names(control$num[control$num])
      new[, nms] <- sweep(new[, nms, drop = FALSE], 2,
                          control$ctr[nms], "-")
      new[, nms] <- sweep(new[, nms, drop = FALSE], 2,
                          control$scl[nms], "/")
    }
    new
  }

  makePreprocWrapper(learner, trainfun, predictfun)
}

## ------------ 4. LEARNERS ------------------------------------------
install.packages("kknn")
library(kknn)
base_learners <- list(
  makeLearner("classif.kknn",          predict.type = "prob"),
  makeLearner("classif.lda",           predict.type = "prob"),
  makeLearner("classif.naiveBayes",    predict.type = "prob"),
  makeLearner("classif.svm",           predict.type = "prob",
              kernel = "radial"),
  makeLearner("classif.rpart",         predict.type = "prob"),
  makeLearner("classif.randomForest",  predict.type = "prob",
              ntree = 500),
  makeLearner("classif.xgboost",       predict.type = "prob",
              nrounds   = 300,
              objective = "multi:softprob")                  # ▲
)
names(base_learners) <- c("knn","lda","nb","svm","rpart","rf","xgb")
wrapped <- lapply(base_learners, makeLimmaWrapper)

## ------------ 5. TASK & RESAMPLING ----------------------------------
task  <- makeClassifTask(id = "liver3",
                         data = data,
                         target = target_var)              # ▲ no “positive”

rdesc <- makeResampleDesc("RepCV", folds = 5, reps = 3, stratify = TRUE)

meas  <- list(multiclass.au1p,        # pairwise AUC (prior-weighted) ▲
              ber,                   # balanced error rate
              mmce,                  # mis-classification error
              kappa)                 # Cohen’s Kappa

## ------------ 6. CROSS-VALIDATION -----------------------------------
results <- lapply(names(wrapped), function(nm) {
  cat("Running", nm, "...\n")
  res <- resample(wrapped[[nm]], task, rdesc,
                  measures = meas, show.info = FALSE)
  data.frame(model  = nm,
             AU1P   = mean(res$measures.test$multiclass.au1p),
             BER    = mean(res$measures.test$ber),
             Error  = mean(res$measures.test$mmce),
             Kappa  = mean(res$measures.test$kappa))
})

cat("\n=======  5×3 cross-validated performance  ========================\n")
print(do.call(rbind, results)[order(-sapply(results, `[[`, "AU1P")), ],
      row.names = FALSE, digits = 3)

```

## Plot the Model comparison
```{r}
## ------------ 7. PLOT PERFORMANCE METRICS ---------------------------
## ------------ 7.  RESAMPLE & COLLECT --------------------------------
summaries <- lapply(names(wrapped), function(nm) {
  res <- resample(wrapped[[nm]], task, rdesc,
                  measures = list(multiclass.au1p, kappa, ber, mmce),
                  show.info = FALSE)
  data.frame(
    model       = nm,
    AU1PMean    = mean(res$measures.test$multiclass.au1p),
    AU1PSD      = sd  (res$measures.test$multiclass.au1p),
    KappaMean   = mean(res$measures.test$kappa),
    KappaSD     = sd  (res$measures.test$kappa),
    BERMean     = mean(res$measures.test$ber),
    BERSD       = sd  (res$measures.test$ber),
    ErrorMean   = mean(res$measures.test$mmce),
    ErrorSD     = sd  (res$measures.test$mmce)
  )
})
res_df <- do.call(rbind, summaries)

## ------------ 8.  TIDY FOR GGPLOT -----------------------------------
library(dplyr)
library(tidyr)

res_long <- res_df %>%
  pivot_longer(-model,
               names_to  = c("Metric", ".value"),
               names_pattern = "(.*)(Mean|SD)")

bar_metrics  <- c("AU1P", "Kappa")      # ▼ bars
bars_long    <- filter(res_long, Metric %in% bar_metrics)

line_metrics <- c("BER", "Error")       # ▼ lines
lines_long   <- filter(res_long, Metric %in% line_metrics)

## ------------ 9.  PLOT ----------------------------------------------
library(ggplot2)
library(colorspace)
library(scales)

## --- colours --------------------------------------------------------
bar_cols <- c(
  "AU1P"  = "#045a8d",
  "Kappa" = "#74a9cf"
)
bar_cols_dark <- darken(bar_cols, 0.25)

line_cols <- c(
  "BER"   = "#ce1256",
  "Error" = "gold3"
)

## --- dodge keeps bars + their error bars grouped --------------------
pos <- position_dodge(width = 0.7)

ggplot() +
  ## BAR metrics ------------------------------------------------------
  geom_col(data = bars_long,
           aes(model, Mean, fill = Metric),
           position = pos, width = 0.65) +
  geom_errorbar(data = bars_long,
                aes(model,
                    ymin = pmax(0, Mean - SD),
                    ymax = pmin(1, Mean + SD),
                    colour = Metric),
                position = pos,
                width = .25, size = .8) +

  ## LINE metrics -----------------------------------------------------
  geom_line(data = lines_long,
            aes(model, Mean, colour = Metric, group = Metric),
            size = 1.2, alpha = 0.6) +
  geom_point(data = lines_long,
             aes(model, Mean, colour = Metric),
             size = 2.2, shape = 16) +
  geom_errorbar(data = lines_long,
                aes(model,
                    ymin = pmax(0, Mean - SD),
                    ymax = pmin(1, Mean + SD),
                    colour = Metric),
                width = .25, size = .8) +

  ## scales -----------------------------------------------------------
  scale_y_continuous(limits = c(0, 1),
                     labels = percent_format(accuracy = 1),
                     name   = "Score (%)") +
  scale_fill_manual(values = bar_cols,   name = NULL) +
  scale_colour_manual(values = c(bar_cols_dark, line_cols), name = NULL) +

  ## theme ------------------------------------------------------------
  theme_minimal(base_size = 12) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  ) +
  scale_x_discrete(
    labels = c(
      knn   = "k-Nearest Neighbors",
      lda   = "Linear Discriminant Analysis (LDA)",
      nb    = "Naïve Bayes",
      svm   = "SVM (RBF)",
      rpart = "Classification And Regression Trees (CART)",
      rf    = "Random Forest",
      xgb   = "XGBoost"
    )
  ) +
  labs(
    title = "5×3 Cross-Validated Performance (mean ± SD)",
    x     = "Model"
  )

## Extract the important features
# pick the wrapper you already built
lrn_xgb <- wrapped[["xgb"]]

# train once on the complete task (no resampling here)
mod_xgb <- mlr::train(lrn_xgb, task)
mod_xgb$factor.levels
# TRUE = keep unwrapping until we hit the learner that owns a real model slot
bst <- getLearnerModel(mod_xgb, more.unwrap = TRUE)
# bst is now an xgboost.Booster

class(bst)

library(xgboost)

# gain = total improvement in loss brought by splits using this feature
imp <- xgb.importance(model = bst)

head(imp, 10)      # top 10 features
install.packages("Ckmeans.1d.dp")   # from CRAN
# (Bioconductor users: choose a CRAN mirror if the prompt appears)

# 2. Load it (xgb.ggplot.importance() checks with requireNamespace)
library(Ckmeans.1d.dp)
xgb.ggplot.importance(imp[1:40], rel_to_first = TRUE)   # bar chart (needs ggplot2)


nclass <- 3
nrounds <- 10

# helper to pick tree indices belonging to class k (0-based)
trees_for_class <- function(k) seq(k, by = nclass, length.out = nrounds)
trees_for_class(0)  # Control
trees_for_class(1)  # Cirrhosis
trees_for_class(2)  # HCC

imp_Control      <- xgb.importance(model = bst, trees = trees_for_class(0))
xgb.ggplot.importance(imp_Control, rel_to_first = TRUE)   # bar chart (needs ggplot2)

imp_Cirrhosis <- xgb.importance(model = bst, trees = trees_for_class(1))
xgb.ggplot.importance(imp_Cirrhosis, rel_to_first = TRUE)

imp_HCC   <- xgb.importance(model = bst, trees = trees_for_class(2))
xgb.ggplot.importance(imp_HCC, rel_to_first = TRUE)

# 1) SHAP values – global & local
install.packages("SHAPforxgboost")
library(SHAPforxgboost)
## 1.1  Pull out the feature table as a *data.frame*
X_df <- as.data.frame(getTaskData(task))

## 1.2  Keep exactly the columns the booster expects
good <- intersect(bst$feature_names, names(X_df))   # sanity check
bad  <- setdiff(bst$feature_names, names(X_df))     # should be character(0)

if (length(bad))
  stop("These ", length(bad), " features are missing in X_df:\n",
       paste(bad, collapse = ", "))

X_df <- X_df[ , good, drop = FALSE]                 # NEVER drop = TRUE here

## 1.3  Convert to plain numeric matrix *with* dim & colnames
X_mat <- data.matrix(X_df)                          # data.matrix keeps colnames 

dim(X_mat)        # e.g. 150   85
head(colnames(X_mat))

library(SHAPforxgboost)
sh <- shap.values(xgb_model = bst, X_train = X_mat)
colnames(X_mat)
# 2) iml::FeatureImp for permutation importance
library(iml)
predictor <- Predictor$new(
  model      = bst,
  data       = X,
  y          = getTaskTargets(task),
  class      = "HCC"  # one class at a time
)
fi <- FeatureImp$new(predictor, loss = "ce")  # cross-entropy
plot(fi)

library(fastshap)
library(mlr)

# 2. Train your final learner on the full dataset


## rerun train() – execution will pause *inside* the xgboost interface
final_mod <- mlr::train(wrapped[["knn"]], task)

# 3. Extract the feature matrix (no target column)
X   <- getTaskData(task, target = FALSE)
X
# 4. Define a prediction wrapper that returns the probability of the positive class
pred_fun <- function(object, newdata) {
  # `object` is the mlr model; we pull out the “prob.diseases” column
  predict(object, newdata = newdata)$data$prob.diseases
}

test <- pred_fun(final_mod, X)
test

# 5. Compute approximate SHAP values
  
library(doParallel)
registerDoParallel(cores = 4)  # use forking with 12 cores
set.seed(2025)

library(future)
library(future.apply)
plan(multisession, workers = availableCores() - 1)  # Windows-friendly

plan(multisession, workers = 4)

system.time({
shap_vals <- fastshap::explain(
  object       = final_mod,     # your trained mlr model
  X            = X,             # background data for conditional sampling
  pred_wrapper = pred_fun,      # must return a numeric vector
  nsim         = 2,  
  adjust = TRUE, parallel = TRUE
             # Monte Carlo reps (tweak higher for stability)
)
})

system.time({
shap_vals <- fastshap::explain(
  object       = final_mod,     # your trained mlr model
  X            = X,             # background data for conditional sampling
  pred_wrapper = pred_fun,      # must return a numeric vector
  nsim         = 2, shap_only = FALSE, 
  adjust = TRUE, parallel = TRUE
             # Monte Carlo reps (tweak higher for stability)
)
})


system.time({
shap_vals_5 <- fastshap::explain(
  object       = final_mod,     # your trained mlr model
  X            = X,             # background data for conditional sampling
  pred_wrapper = pred_fun,      # must return a numeric vector
  newdata = X[5,],
  nsim         = 5,  
  adjust = TRUE, parallel = TRUE
             # Monte Carlo reps (tweak higher for stability)
)
})

saveRDS(shap_vals, file = "/Users/macos/Desktop/AFLATOXINS/shap_vals.rds" )
baseline <- attr(shap_vals, "baseline")

install.packages("shapviz")
library(shapviz)

# Local explanations
## Change the compound name
shv <- shapviz(shap_vals_5, X = X[5,], baseline = baseline)

id2name <- rowData(se_50)$Description                # or whatever your column is called
names(id2name) <- paste0("meta_", rownames(se_50)) # prepend “meta_” to match shap names
view(rowData(se_50))
## ------------------------------------------------------------------
## 2.  Identify the features present in the shapviz object
## ------------------------------------------------------------------
old_names <- colnames(X)                # SHAP matrix column names
hits      <- intersect(old_names, names(id2name)) # those we can translate
length(hits)
## --- update SHAP matrix dimnames -------------------------------------
col_idx <- match(hits, old_names)
ok      <- which(!is.na(col_idx))
colnames(X)[col_idx[ok]] <- id2name[hits[ok]]

## --- update shapviz 'feature_names' attribute ------------------------
feat_idx <- match(hits, attr(shv, "feature_names"))
ok2      <- which(!is.na(feat_idx))
attr(shv, "feature_names")[feat_idx[ok2]] <- id2name[hits[ok2]]


## 2. give it the same names we assigned to S

X_new <- X                      # original data you used to build shv
colnames(X_new)[match(hits, colnames(X_new))[ok]] <- id2name[hits[ok]]
dim(X_new)
dim(shv)

colnames(X_new[5,])
colnames(shv)
colnames(shv$S) 
identical(colnames(X_new[5,]), colnames(shv))
identical(colnames(shv$S), colnames(X_new)) 

sv_waterfall(shv, fill_colors = c("#a52c60", "#f7d13d"))


sv_force(shv, fill_colors = c("#a52c60", "#f7d13d"))
# 6. Aggregate into mean absolute contributions
mean_abs_shap <- colMeans(abs(shap_vals))
importance <- data.frame(
  Variable      = names(mean_abs_shap),
  MeanAbsSHAP   = mean_abs_shap
)
importance <- importance[order(importance$MeanAbsSHAP, decreasing = TRUE), ]

# 7. View the ranked list
print(importance)

# Global explanations
baseline <- attr(shap_vals, "baseline")
tshv.global <- shapviz(shap_vals, X = X, baseline = baseline)
sv_importance(tshv.global, max_display = 30L, fill = "#e08214")  

sv_dependence(tshv.global, v = "bmi_new", alpha = 0.3)

