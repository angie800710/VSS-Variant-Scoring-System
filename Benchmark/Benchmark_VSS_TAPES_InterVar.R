library(dplyr)
library(purrr)
library(tibble)
library(DescTools)

setwd("C:/Users/User1/OneDrive/桌面/VSS2_R/Revision")

### Import and merge data ###
merge <- read.csv("Merge_vss_tapes_intervar.csv", header = T)


### ROC compare: TAPES vs. VSS
roc_data <- merge %>%
  mutate(
    Truth_binary = case_when(
      Assertion %in% c(
        "Pathogenic",
        "Likely Pathogenic"
      ) ~ 1,
      
      Assertion %in% c(
        "Benign",
        "Likely Benign"
      ) ~ 0,
      
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Truth_binary),
    !is.na(VSS2_score),
    !is.na(Tapes_prob)
  )

# Check sample size
nrow(roc_data) #8746
table(roc_data$Truth_binary)


### ============================================================
### 2. Calculate ROC curves
### ============================================================

roc_vss <- roc(
  response = roc_data$Truth_binary,
  predictor = roc_data$VSS2_score,
  levels = c(0, 1),
  direction = "<",
  ci = TRUE
)

roc_tapes <- roc(
  response = roc_data$Truth_binary,
  predictor = roc_data$Tapes_prob,
  levels = c(0, 1),
  direction = "<",
  ci = TRUE
)


### ============================================================
### 3. AUROC and 95% CI
### ============================================================

auc(roc_vss)
ci.auc(roc_vss)

auc(roc_tapes)
ci.auc(roc_tapes)

roc.test(
  roc_vss,
  roc_tapes,
  method = "delong",
  paired = TRUE
)

auc_vss   <- as.numeric(auc(roc_vss)) #0.99
auc_tapes <- as.numeric(auc(roc_tapes)) #0.96

coords(roc_vss, "best", ret=c("threshold","specificity", "sensitivity"), 
       best.method=c("youden", "closest.topleft"))
coords(roc_tapes, "best", ret=c("threshold","specificity", "sensitivity"), 
       best.method=c("youden", "closest.topleft"))


### ============================================================
### 4. Paired DeLong test
### ============================================================

roc_test <- roc.test(
  roc_vss,
  roc_tapes,
  paired = TRUE,
  method = "delong"
)

roc_test


### ============================================================
### 5. ROC comparison plot
### ============================================================

# ROC plot for pathogenic prediction
### ============================================================
### Extract ROC coordinates
### ============================================================

roc_vss_df <- pROC::coords(
  roc_vss,
  x = "all",
  ret = c("specificity", "sensitivity"),
  transpose = FALSE
) %>%
  as.data.frame() %>%
  mutate(
    FPR = 1 - specificity,
    TPR = sensitivity,
    Tool = "VSS"
  )

roc_tapes_df <- pROC::coords(
  roc_tapes,
  x = "all",
  ret = c("specificity", "sensitivity"),
  transpose = FALSE
) %>%
  as.data.frame() %>%
  mutate(
    FPR = 1 - specificity,
    TPR = sensitivity,
    Tool = "TAPES"
  )

roc_df <- bind_rows(
  roc_vss_df,
  roc_tapes_df
)


### ============================================================
### Plot ROC curves
### ============================================================
roc_df$Tool <- factor(
  roc_df$Tool,
  levels = c("VSS", "TAPES")
)


p_roc_final <- ggplot(
  roc_df,
  aes(
    x = FPR,
    y = TPR,
    color = Tool
  )
) +
  geom_line(linewidth = 1) +
  
  # No-discrimination reference line
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "VSS" = "red",
      "TAPES" = "blue"
    ),
    labels = c(
      "VSS" = paste0(
        "VSS (AUROC = ",
        sprintf("%.3f", auc_vss),
        ")"
      ),
      "TAPES" = paste0(
        "TAPES (AUROC = ",
        sprintf("%.3f", auc_tapes),
        ")"
      )
    )
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  
  coord_equal() +
  
  labs(
    x = "1 − Specificity",
    y = "Sensitivity"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.title = element_blank(),
    legend.position = c(0.72, 0.18)
  )

p_roc_final
  
  
  ### ============================================================
  ### Plot ROC curves
  ### ============================================================
  
  p_roc_final <- ggplot(
    roc_df,
    aes(
      x = FPR,
      y = TPR,
      color = Tool
    )
  ) +
    geom_line(linewidth = 1) +
    
    # No-discrimination reference line
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
    scale_color_manual(
      values = c(
        "VSS" = "red",
        "TAPES" = "blue"
      ),
      labels = c(
        "VSS" = paste0(
          "VSS (AUROC = ",
          sprintf("%.3f", auc_vss),
          ")"
        ),
        "TAPES" = paste0(
          "TAPES (AUROC = ",
          sprintf("%.3f", auc_tapes),
          ")"
        )
      )
    ) +
    
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    
    coord_equal() +
    
    labs(
      x = "1 − Specificity",
      y = "Sensitivity"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      legend.title = element_blank(),
      legend.position = c(0.72, 0.18)
    )
  
  p_roc_final
  

  ### ============================================================
  ### 1. Calculate PR curves
  ### Positive = P/LP
  ### Negative = B/LB
  ### ============================================================
  
  pr_vss <- pr.curve(
    scores.class0 = roc_data$VSS2_score[
      roc_data$Truth_binary == 1
    ],
    scores.class1 = roc_data$VSS2_score[
      roc_data$Truth_binary == 0
    ],
    curve = TRUE
  )
  
  pr_tapes <- pr.curve(
    scores.class0 = roc_data$Tapes_prob[
      roc_data$Truth_binary == 1
    ],
    scores.class1 = roc_data$Tapes_prob[
      roc_data$Truth_binary == 0
    ],
    curve = TRUE
  )
  
  
  ### ============================================================
  ### 2. AUPRC
  ### ============================================================
  
  auprc_vss <- pr_vss$auc.integral
  auprc_tapes <- pr_tapes$auc.integral
  
  auprc_vss
  auprc_tapes
  
  ### ============================================================
  ### 檢定 Paired bootstrap test for AUPRC
  ### ============================================================
  
  set.seed(1234)
  
  B <- 5000
  n <- nrow(roc_data)
  
  boot_auprc_diff <- replicate(B, {
    
    # Paired bootstrap: resample the same variants for both tools
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    boot_data <- roc_data[idx, ]
    
    # VSS
    pr_vss_boot <- pr.curve(
      scores.class0 = boot_data$VSS2_score[
        boot_data$Truth_binary == 1
      ],
      scores.class1 = boot_data$VSS2_score[
        boot_data$Truth_binary == 0
      ],
      curve = FALSE
    )
    
    # TAPES
    pr_tapes_boot <- pr.curve(
      scores.class0 = boot_data$Tapes_prob[
        boot_data$Truth_binary == 1
      ],
      scores.class1 = boot_data$Tapes_prob[
        boot_data$Truth_binary == 0
      ],
      curve = FALSE
    )
    
    # Difference: VSS - TAPES
    pr_vss_boot$auc.integral -
      pr_tapes_boot$auc.integral
  })
  
  auprc_diff <- auprc_vss - auprc_tapes
  auprc_diff
  
  auprc_diff_ci <- quantile(
    boot_auprc_diff,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  auprc_diff_ci
  
  p_auprc <- 2 * min(
    mean(boot_auprc_diff <= 0, na.rm = TRUE),
    mean(boot_auprc_diff >= 0, na.rm = TRUE)
  )
  p_auprc
  
  cat(
    "VSS AUPRC =", round(auprc_vss, 4), "\n",
    "TAPES AUPRC =", round(auprc_tapes, 4), "\n",
    "Difference (VSS - TAPES) =", round(auprc_diff, 4), "\n",
    "95% bootstrap CI =",
    round(auprc_diff_ci[1], 4), "to",
    round(auprc_diff_ci[2], 4), "\n",
    "Paired bootstrap P =",
    format.pval(p_auprc, digits = 3)
  )
  
  
  
  ### ============================================================
  ### 3. Extract PR curve coordinates
  ### ============================================================
  
  pr_vss_df <- as.data.frame(pr_vss$curve)
  colnames(pr_vss_df) <- c(
    "Recall",
    "Precision",
    "Threshold"
  )
  
  pr_vss_df <- pr_vss_df %>%
    mutate(
      Tool = "VSS"
    )
  
  
  pr_tapes_df <- as.data.frame(pr_tapes$curve)
  colnames(pr_tapes_df) <- c(
    "Recall",
    "Precision",
    "Threshold"
  )
  
  pr_tapes_df <- pr_tapes_df %>%
    mutate(
      Tool = "TAPES"
    )
  
  
  ### Combine
  pr_df <- bind_rows(
    pr_vss_df,
    pr_tapes_df
  ) %>%
    mutate(
      Tool = factor(
        Tool,
        levels = c("VSS", "TAPES")
      )
    )
  
  
  ### ============================================================
  ### 4. Positive prevalence
  ### No-skill baseline for PR curve
  ### ============================================================
  
  positive_prevalence <- mean(
    roc_data$Truth_binary == 1,
    na.rm = TRUE
  )
  
  positive_prevalence
  
  
  ### ============================================================
  ### 5. Plot PR curve
  ### ============================================================
  
  p_pr_final <- ggplot(
    pr_df,
    aes(
      x = Recall,
      y = Precision,
      color = Tool
    )
  ) +
    
    geom_line(
      linewidth = 1
    ) +
    
    # No-skill baseline
    geom_hline(
      yintercept = positive_prevalence,
      linetype = "dashed",
      linewidth = 0.7
    ) +
    
    scale_color_manual(
      values = c(
        "VSS" = "red",
        "TAPES" = "blue"
      ),
      breaks = c(
        "VSS",
        "TAPES"
      ),
      labels = c(
        "VSS" = paste0(
          "VSS (AUPRC = ",
          sprintf("%.3f", auprc_vss),
          ")"
        ),
        "TAPES" = paste0(
          "TAPES (AUPRC = ",
          sprintf("%.3f", auprc_tapes),
          ")"
        )
      )
    ) +
    
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    
    coord_equal() +
    
    labs(
      x = "Recall",
      y = "Precision"
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      legend.title = element_blank(),
      legend.position = c(0.72, 0.18)
    )
  
  p_pr_final
  
  
### ============================================================
### 3-class benchmark
### Reference: ClinGen
### Classes: P/LP, VUS, B/LB
### ============================================================

benchmark_3class <- merge %>%
  mutate(
    
    ## ClinGen ground truth
    Truth = case_when(
      Assertion %in% c("Pathogenic", "Likely Pathogenic") ~ "P/LP",
      
      Assertion == "Uncertain Significance" ~ "VUS",
      
      Assertion %in% c("Benign", "Likely Benign") ~ "B/LB",
      
      TRUE ~ NA_character_
    ),
    
    
    ## VSS prediction
    ## P/LP >= 6
    ## VUS = 2–5
    ## B/LB <= 1
    VSS_pred = case_when(
      is.na(VSS2_score) ~ NA_character_,
      VSS2_score >= 6 ~ "P/LP",
      VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
      VSS2_score <= 1 ~ "B/LB"
    ),
    
    
    ## TAPES prediction
    TAPES_pred = case_when(
      is.na(Tapes_prediction) ~ NA_character_,
      
      Tapes_prediction %in%
        c("Pathogenic", "Likely Pathogenic") ~ "P/LP",
      
      Tapes_prediction == "VUS" ~ "VUS",
      
      Tapes_prediction %in%
        c("Benign", "Likely Benign") ~ "B/LB",
      
      TRUE ~ NA_character_
    ),
    
    
    ## InterVar prediction
    InterVar_pred = case_when(
      is.na(InterVar_classification) ~ NA_character_,
      
      InterVar_classification %in%
        c("Pathogenic", "Likely pathogenic") ~ "P/LP",
      
      InterVar_classification ==
        "Uncertain significance" ~ "VUS",
      
      InterVar_classification %in%
        c("Benign", "Likely benign") ~ "B/LB",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Retain only variants with ClinGen classification
  filter(!is.na(Truth))


### ============================================================
### Restrict to common variants classified by all three tools
### ============================================================

benchmark_3class_common <- benchmark_3class %>%
  filter(
    !is.na(VSS_pred),
    !is.na(TAPES_pred),
    !is.na(InterVar_pred)
  )

# Check sample size and ClinGen class distribution
nrow(benchmark_3class_common) #11414
table(benchmark_3class_common$Truth) #P/VUS/B:5770/2839/2805


### ============================================================
### Function for multiclass performance metrics
### One-vs-rest metrics + macro-average
### ============================================================

calc_3class_metrics <- function(data, pred_col, tool_name) {
  
  classes <- c("P/LP", "VUS", "B/LB")
  
  # Overall multiclass accuracy
  overall_accuracy <- mean(
    data$Truth == data[[pred_col]]
  )
  
  # One-vs-rest metrics for each class
  per_class <- map_dfr(classes, function(cls) {
    
    truth_positive <- data$Truth == cls
    pred_positive  <- data[[pred_col]] == cls
    
    TP <- sum(truth_positive & pred_positive)
    TN <- sum(!truth_positive & !pred_positive)
    FP <- sum(!truth_positive & pred_positive)
    FN <- sum(truth_positive & !pred_positive)
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    precision   <- TP / (TP + FP)
    recall      <- sensitivity
    ppv         <- precision
    npv         <- TN / (TN + FN)
    
    f1_score <- ifelse(
      is.na(precision) |
        is.na(recall) |
        (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)
    )
    
    tibble(
      Tool = tool_name,
      Class = cls,
      N = nrow(data),
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      Recall = recall,
      PPV = ppv,
      NPV = npv,
      F1_score = f1_score,
      Accuracy = overall_accuracy
    )
  })
  
  
  # Macro-average across the three classes
  macro_average <- per_class %>%
    summarise(
      Tool = tool_name,
      Class = "Macro-average",
      N = first(N),
      TP = NA_integer_,
      TN = NA_integer_,
      FP = NA_integer_,
      FN = NA_integer_,
      Sensitivity = mean(Sensitivity, na.rm = TRUE),
      Specificity = mean(Specificity, na.rm = TRUE),
      Precision = mean(Precision, na.rm = TRUE),
      Recall = mean(Recall, na.rm = TRUE),
      PPV = mean(PPV, na.rm = TRUE),
      NPV = mean(NPV, na.rm = TRUE),
      F1_score = mean(F1_score, na.rm = TRUE),
      Accuracy = first(Accuracy)
    )
  
  bind_rows(per_class, macro_average)
}


### ============================================================
### Calculate metrics for VSS, TAPES and InterVar
### ============================================================

metrics_3class_common <- bind_rows(
  
  calc_3class_metrics(
    benchmark_3class_common,
    "VSS_pred",
    "VSS"
  ),
  
  calc_3class_metrics(
    benchmark_3class_common,
    "TAPES_pred",
    "TAPES"
  ),
  
  calc_3class_metrics(
    benchmark_3class_common,
    "InterVar_pred",
    "InterVar"
  )
)


### ============================================================
### Final performance comparison table
### ============================================================

performance_compare_3class <- metrics_3class_common %>%
  filter(Class == "Macro-average") %>%
  transmute(
    Tool,
    N,
    Accuracy = round(Accuracy * 100, 2),
    Macro_Precision = round(Precision * 100, 2),
    Macro_Recall = round(Recall * 100, 2),
    Macro_F1_score = round(F1_score * 100, 2)
  )

performance_compare_3class


### ============================================================
### 2-class benchmark
### Reference: ClinGen
### Classes: P/LP vs B/LB
### VUS excluded
### VSS cutoff: >= 6 = P/LP
### ============================================================

benchmark_2class <- merge %>%
  mutate(
    
    ## ClinGen ground truth
    Truth = case_when(
      Assertion %in% c(
        "Pathogenic",
        "Likely Pathogenic"
      ) ~ "P/LP",
      
      Assertion %in% c(
        "Benign",
        "Likely Benign"
      ) ~ "B/LB",
      
      TRUE ~ NA_character_
    ),
    
    
    ## VSS prediction
    ## P/LP >= 6
    ## B/LB < 6
    VSS_pred = case_when(
      is.na(VSS2_score) ~ NA_character_,
      VSS2_score >= 6 ~ "P/LP",
      VSS2_score < 6 ~ "B/LB"
    ),
    
    
    ## TAPES prediction
    TAPES_pred = case_when(
      is.na(Tapes_prediction) ~ NA_character_,
      
      Tapes_prediction %in%
        c("Pathogenic", "Likely Pathogenic") ~ "P/LP",
      
      Tapes_prediction %in%
        c("Benign", "Likely Benign") ~ "B/LB",
      
      ## TAPES VUS excluded
      TRUE ~ NA_character_
    ),
    
    
    ## InterVar prediction
    InterVar_pred = case_when(
      is.na(InterVar_classification) ~ NA_character_,
      
      InterVar_classification %in%
        c("Pathogenic", "Likely pathogenic") ~ "P/LP",
      
      InterVar_classification %in%
        c("Benign", "Likely benign") ~ "B/LB",
      
      ## InterVar VUS excluded
      TRUE ~ NA_character_
    )
  ) %>%
  
  ## Remove ClinGen VUS and other non-binary assertions
  filter(!is.na(Truth))


### ============================================================
### Restrict to variants classified by all three tools
### ============================================================

benchmark_2class_common <- benchmark_2class %>%
  filter(
    !is.na(VSS_pred),
    !is.na(TAPES_pred),
    !is.na(InterVar_pred)
  )

# Sample size
nrow(benchmark_2class_common)

# ClinGen class distribution
table(benchmark_2class_common$Truth)

# Predictions by each tool
table(benchmark_2class_common$VSS_pred)
table(benchmark_2class_common$TAPES_pred)
table(benchmark_2class_common$InterVar_pred)


### Function for binary performance metrics

calc_2class_metrics <- function(data, pred_col, tool_name) {
  
  classes <- c("P/LP", "B/LB")
  
  # Overall multiclass/binary accuracy
  overall_accuracy <- mean(
    data$Truth == data[[pred_col]]
  )
  
  
  ### One-vs-rest metrics for each class
  per_class <- map_dfr(classes, function(cls) {
    
    truth_positive <- data$Truth == cls
    pred_positive  <- data[[pred_col]] == cls
    
    TP <- sum(truth_positive & pred_positive)
    TN <- sum(!truth_positive & !pred_positive)
    FP <- sum(!truth_positive & pred_positive)
    FN <- sum(truth_positive & !pred_positive)
    
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    
    precision <- TP / (TP + FP)
    recall <- sensitivity
    
    ppv <- precision
    npv <- TN / (TN + FN)
    
    f1_score <- ifelse(
      is.na(precision) |
        is.na(recall) |
        (precision + recall) == 0,
      NA_real_,
      2 * precision * recall / (precision + recall)
    )
    
    tibble(
      Tool = tool_name,
      Class = cls,
      N = nrow(data),
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      Recall = recall,
      PPV = ppv,
      NPV = npv,
      F1_score = f1_score,
      Accuracy = overall_accuracy
    )
  })
  
  
  ### Macro-average across two classes
  macro_average <- per_class %>%
    summarise(
      Tool = tool_name,
      Class = "Macro-average",
      N = first(N),
      TP = NA_integer_,
      TN = NA_integer_,
      FP = NA_integer_,
      FN = NA_integer_,
      Sensitivity = mean(Sensitivity, na.rm = TRUE),
      Specificity = mean(Specificity, na.rm = TRUE),
      Precision = mean(Precision, na.rm = TRUE),
      Recall = mean(Recall, na.rm = TRUE),
      PPV = mean(PPV, na.rm = TRUE),
      NPV = mean(NPV, na.rm = TRUE),
      F1_score = mean(F1_score, na.rm = TRUE),
      Accuracy = first(Accuracy)
    )
  
  bind_rows(
    per_class,
    macro_average
  )
}

metrics_2class_common <- bind_rows(
  
  calc_2class_metrics(
    benchmark_2class_common,
    "VSS_pred",
    "VSS"
  ),
  
  calc_2class_metrics(
    benchmark_2class_common,
    "TAPES_pred",
    "TAPES"
  ),
  
  calc_2class_metrics(
    benchmark_2class_common,
    "InterVar_pred",
    "InterVar"
  )
)

metrics_2class_common

performance_compare_2class <- metrics_2class_common %>%
  filter(Class == "Macro-average") %>%
  transmute(
    Tool,
    N,
    Accuracy = round(Accuracy * 100, 2),
    Macro_Precision = round(Precision * 100, 2),
    Macro_Recall = round(Recall * 100, 2),
    Macro_F1_score = round(F1_score * 100, 2)
  )

performance_compare_2class

### ============================================================
### 5-class benchmark
### Reference: ClinGen
### Classes: P, LP, VUS, LB, B
### ============================================================

benchmark_5class <- merge %>%
  mutate(
    
    ## ClinGen ground truth
    Truth = case_when(
      Assertion == "Pathogenic" ~ "P",
      
      Assertion == "Likely Pathogenic" ~ "LP",
      
      Assertion == "Uncertain Significance" ~ "VUS",
      
      Assertion == "Likely Benign" ~ "LB",
      
      Assertion == "Benign" ~ "B",
      
      TRUE ~ NA_character_
    ),
    
    
    ## VSS prediction
    ## P >= 10
    ## LP = 6–9
    ## VUS = 2–5
    ## LB = 0–1
    ## B < 0
    VSS_pred = case_when(
      is.na(VSS2_score) ~ NA_character_,
      
      VSS2_score >= 10 ~ "P",
      
      VSS2_score >= 6 & VSS2_score <= 9 ~ "LP",
      
      VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
      
      VSS2_score >= 0 & VSS2_score <= 1 ~ "LB",
      
      VSS2_score < 0 ~ "B"
    ),
    
    
    ## TAPES prediction
    TAPES_pred = case_when(
      is.na(Tapes_prediction) ~ NA_character_,
      
      Tapes_prediction == "Pathogenic" ~ "P",
      
      Tapes_prediction == "Likely Pathogenic" ~ "LP",
      
      Tapes_prediction == "VUS" ~ "VUS",
      
      Tapes_prediction == "Likely Benign" ~ "LB",
      
      Tapes_prediction == "Benign" ~ "B",
      
      TRUE ~ NA_character_
    ),
    
    
    ## InterVar prediction
    InterVar_pred = case_when(
      is.na(InterVar_classification) ~ NA_character_,
      
      InterVar_classification == "Pathogenic" ~ "P",
      
      InterVar_classification == "Likely pathogenic" ~ "LP",
      
      InterVar_classification ==
        "Uncertain significance" ~ "VUS",
      
      InterVar_classification == "Likely benign" ~ "LB",
      
      InterVar_classification == "Benign" ~ "B",
      
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Retain only variants with ClinGen classification
  filter(!is.na(Truth))


### ============================================================
### Restrict to common variants classified by all three tools
### ============================================================

benchmark_5class_common <- benchmark_5class %>%
  filter(
    !is.na(VSS_pred),
    !is.na(TAPES_pred),
    !is.na(InterVar_pred)
  )

# Check sample size and ClinGen class distribution
nrow(benchmark_5class_common)

table(
  benchmark_5class_common$Truth
)


### ============================================================
### Function for 5-class performance metrics
### One-vs-rest metrics + macro-average
### ============================================================

calc_5class_metrics <- function(data, pred_col, tool_name) {
  
  classes <- c(
    "P",
    "LP",
    "VUS",
    "LB",
    "B"
  )
  
  # Overall multiclass accuracy
  overall_accuracy <- mean(
    data$Truth == data[[pred_col]]
  )
  
  
  # One-vs-rest metrics for each class
  per_class <- map_dfr(
    classes,
    function(cls) {
      
      truth_positive <- data$Truth == cls
      pred_positive  <- data[[pred_col]] == cls
      
      TP <- sum(
        truth_positive & pred_positive
      )
      
      TN <- sum(
        !truth_positive & !pred_positive
      )
      
      FP <- sum(
        !truth_positive & pred_positive
      )
      
      FN <- sum(
        truth_positive & !pred_positive
      )
      
      
      sensitivity <- TP / (TP + FN)
      specificity <- TN / (TN + FP)
      
      precision <- TP / (TP + FP)
      recall <- sensitivity
      
      ppv <- precision
      npv <- TN / (TN + FN)
      
      f1_score <- ifelse(
        is.na(precision) |
          is.na(recall) |
          (precision + recall) == 0,
        NA_real_,
        2 * precision * recall /
          (precision + recall)
      )
      
      
      tibble(
        Tool = tool_name,
        Class = cls,
        N = nrow(data),
        TP = TP,
        TN = TN,
        FP = FP,
        FN = FN,
        Sensitivity = sensitivity,
        Specificity = specificity,
        Precision = precision,
        Recall = recall,
        PPV = ppv,
        NPV = npv,
        F1_score = f1_score,
        Accuracy = overall_accuracy
      )
    }
  )
  
  
  ### Macro-average across five classes
  macro_average <- per_class %>%
    summarise(
      Tool = tool_name,
      Class = "Macro-average",
      N = first(N),
      TP = NA_integer_,
      TN = NA_integer_,
      FP = NA_integer_,
      FN = NA_integer_,
      Sensitivity = mean(
        Sensitivity,
        na.rm = TRUE
      ),
      Specificity = mean(
        Specificity,
        na.rm = TRUE
      ),
      Precision = mean(
        Precision,
        na.rm = TRUE
      ),
      Recall = mean(
        Recall,
        na.rm = TRUE
      ),
      PPV = mean(
        PPV,
        na.rm = TRUE
      ),
      NPV = mean(
        NPV,
        na.rm = TRUE
      ),
      F1_score = mean(
        F1_score,
        na.rm = TRUE
      ),
      Accuracy = first(Accuracy)
    )
  
  
  bind_rows(
    per_class,
    macro_average
  )
}


### ============================================================
### Calculate metrics for VSS, TAPES and InterVar
### ============================================================

metrics_5class_common <- bind_rows(
  
  calc_5class_metrics(
    benchmark_5class_common,
    "VSS_pred",
    "VSS"
  ),
  
  calc_5class_metrics(
    benchmark_5class_common,
    "TAPES_pred",
    "TAPES"
  ),
  
  calc_5class_metrics(
    benchmark_5class_common,
    "InterVar_pred",
    "InterVar"
  )
)

metrics_5class_common


### ============================================================
### Final 5-class performance comparison table
### ============================================================

performance_compare_5class <- metrics_5class_common %>%
  filter(
    Class == "Macro-average"
  ) %>%
  transmute(
    Tool,
    N,
    Accuracy = round(
      Accuracy * 100,
      2
    ),
    Macro_Precision = round(
      Precision * 100,
      2
    ),
    Macro_Recall = round(
      Recall * 100,
      2
    ),
    Macro_F1_score = round(
      F1_score * 100,
      2
    )
  )

performance_compare_5class
