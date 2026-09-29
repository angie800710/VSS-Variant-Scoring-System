library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(furrr)
library(purrr)
library(stringr)
library(reshape2)
library(data.table)
library(DET)
library(readr)
library(pROC)
library(plotROC)
library(PRROC)
library(scales)
library(purrr)
library(tibble)
library(caret)

setwd("C:/Users/User1/OneDrive/桌面/VSS2_R/Revision")

### Import and merge data ###
VSS <- read.csv("VSS2_ClinGen.csv", header = T)
Tapes <- read.csv("Tapes_ClinGen.csv", header = T)
Intervar <- read.csv("InterVar_ClinGen.csv", header = T)

VSS %>%
  count(Index) %>%
  filter(n > 1)

Tapes %>%
  count(Index) %>%
  filter(n > 1)

Intervar %>%
  count(Index) %>%
  filter(n > 1)

Tapes <- Tapes %>%
  distinct(Index, .keep_all = TRUE)

Intervar <- Intervar %>%
  distinct(Index, .keep_all = TRUE)

merge <- VSS %>%
  inner_join(Tapes, by = "Index") %>%
  inner_join(Intervar, by = "Index")

write.csv(merge, "Merge_vss_tapes_intervar.csv", row.names=F)

### ClinGen vs. VSS ###
colnames(VSS) #12619筆

# Remove VUS and create binary ground truth
VSS_binary <- VSS %>%
  filter(Assertion %in% c("Pathogenic", "Likely Pathogenic", "Benign", "Likely Benign")) %>%
  mutate(
    Truth = if_else(Assertion %in% c("Pathogenic", "Likely Pathogenic"), 1, 0)
  ) #8852筆

##ROC curve + AUC##
roc_obj <- roc(
  response = VSS_binary$Truth,
  predictor = VSS_binary$VSS2_score,
  levels = c(0, 1),
  direction = "<"
)

auc(roc_obj)

plot(
  roc_obj,
  print.auc = TRUE,
  legacy.axes = TRUE,
  main = "ROC Curve of VSS"
)

# Optimal threshold
coords(roc_obj, x = 2, input="threshold",
       ret=c("sensitivity","specificity"))
coords(roc_obj, x = 3, input="threshold",
       ret=c("sensitivity","specificity"))
coords(roc_obj, x = 4, input="threshold",
       ret=c("sensitivity","specificity"))
coords(roc_obj, x = 5, input="threshold",
       ret=c("sensitivity","specificity"))
coords(roc_obj, x = 6, input="threshold",
       ret=c("sensitivity","specificity"))

coords(roc_obj, "best", ret=c("threshold","specificity", "sensitivity"), 
       best.method=c("youden", "closest.topleft"))

## PRROC ##
## Separate positive and negative scores ##

# P / LP
scores_pos <- VSS_binary$VSS2_score[VSS_binary$Truth == 1]

# B / LB
scores_neg <- VSS_binary$VSS2_score[VSS_binary$Truth == 0]


## AUROC ##

roc_prroc <- roc.curve(
  scores.class0 = scores_pos,
  scores.class1 = scores_neg,
  curve = TRUE
)

roc_prroc$auc

plot(
  roc_prroc,
  main = paste0(
    "ROC Curve of VSS (AUROC = ",
    round(roc_prroc$auc, 3),
    ")"
  )
)


## AUPRC ##

pr_prroc <- pr.curve(
  scores.class0 = scores_pos,
  scores.class1 = scores_neg,
  curve = TRUE
)

pr_prroc$auc.integral

plot(
  pr_prroc,
  main = paste0(
    "PR Curve of VSS (AUPRC = ",
    round(pr_prroc$auc.integral, 3),
    ")"
  )
)

### Benign ROC ###
# B/LB = 1 (positive)
# P/LP = 0 (negative)
VSS_binary$Truth_B <- ifelse(VSS_binary$Truth == 0, 1, 0)

roc_benign <- roc(
  response = VSS_binary$Truth_B,
  predictor = VSS_binary$VSS2_score,
  levels = c(0, 1),
  direction = ">"
)




## Metrics ##
# Ground truth
true_labels <- VSS_binary$Truth

# VSS possible cutoffs
cutoffs <- seq(-2, 13, 1)


## Calculate metrics at each cutoff ##

calculate_metrics <- function(cutoff, pred_scores, true_labels) {
  
  # VSS score >= cutoff = Pathogenic
  predicted <- ifelse(pred_scores >= cutoff, 1, 0)
  
  # Confusion matrix
  TP <- sum(predicted == 1 & true_labels == 1)
  TN <- sum(predicted == 0 & true_labels == 0)
  FP <- sum(predicted == 1 & true_labels == 0)
  FN <- sum(predicted == 0 & true_labels == 1)
  
  # Metrics
  Sensitivity <- TP / (TP + FN)
  Specificity <- TN / (TN + FP)
  
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  Precision <- PPV
  Recall <- Sensitivity
  
  F1_score <- 2 * Precision * Recall / (Precision + Recall)
  
  # Youden's J / Youden's index
  Youden_J <- Sensitivity + Specificity - 1
  
  return(c(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    Sensitivity = Sensitivity,
    Specificity = Specificity,
    PPV = PPV,
    NPV = NPV,
    Accuracy = Accuracy,
    Precision = Precision,
    Recall = Recall,
    F1_score = F1_score,
    Youden_J = Youden_J
  ))
}


## Run all cutoffs ##

metrics <- t(
  sapply(
    cutoffs,
    calculate_metrics,
    pred_scores = VSS_binary$VSS2_score,
    true_labels = VSS_binary$Truth
  )
)

metrics_df <- as.data.frame(metrics)

metrics_df$Cutoff <- cutoffs

# Put cutoff first
metrics_df <- metrics_df %>%
  select(
    Cutoff,
    TP, TN, FP, FN,
    Sensitivity, Specificity,
    PPV, NPV,
    Youden_J,
    Accuracy,
    Precision, Recall,
    F1_score
  )

metrics_df
write.csv(metrics_df, "VSS_ClinGen_metrics.csv", row.names=F)

## Plot C. SEN/SPE crossover at 5 ##
plot_C <- metrics_df %>%
  select(Cutoff, Sensitivity, Specificity) %>%
  pivot_longer(
    cols = c(Sensitivity, Specificity),
    names_to = "Metric",
    values_to = "Value"
  )

ggplot(
  plot_C,
  aes(x = Cutoff, y = Value, group = Metric, linetype = Metric)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = 5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = seq(-2, 13, 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    x = "VSS score cutoff",
    y = "Performance",
    linetype = NULL
  ) +
  theme_classic()

## Plot D. Histogram分布圖 ##
## Prepare data for Figure 1D ##

VSS_plotD <- VSS %>%
  filter(
    Assertion %in% c(
      "Pathogenic",
      "Likely Pathogenic",
      "Benign",
      "Likely Benign"
    )
  ) %>%
  mutate(
    ClinGen_class = case_when(
      Assertion %in% c("Pathogenic", "Likely Pathogenic") ~ "P/LP",
      Assertion %in% c("Benign", "Likely Benign") ~ "B/LB"
    )
  )


## Calculate score distribution within each category ##

score_dist <- VSS_plotD %>%
  count(ClinGen_class, VSS2_score) %>%
  group_by(ClinGen_class) %>%
  mutate(
    Proportion = n / sum(n)
  ) %>%
  ungroup()


## Complete missing score × category combinations ##
## This keeps all bars the same width even when count = 0 ##

score_dist_complete <- score_dist %>%
  complete(
    VSS2_score = seq(-2, 13, 1),
    ClinGen_class = c("P/LP", "B/LB"),
    fill = list(
      n = 0,
      Proportion = 0
    )
  ) %>%
  mutate(
    ClinGen_class = factor(
      ClinGen_class,
      levels = c("P/LP", "B/LB")
    )
  )

write.csv(score_dist_complete, "VSS_score_prop_ClinGen.csv", row.names=F)

## Figure 1D ##

plot_D <- ggplot(
  score_dist_complete,
  aes(
    x = VSS2_score,
    y = Proportion,
    fill = ClinGen_class
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.85,
      preserve = "single"
    ),
    width = 0.8
  ) +
  
  # Youden's index optimal threshold
  geom_vline(
    xintercept = 5.5,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "P/LP" = "#F8766D",
      "B/LB" = "#00BFC4"
    ),
    breaks = c("P/LP", "B/LB")
  ) +
  
  scale_x_continuous(
    breaks = seq(-2, 13, 1)
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = "VSS score",
    y = "Proportion of variants",
    fill = NULL
  ) +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.position = "right",
    legend.text = element_text(size = 13)
  )

plot_D
##########################Benchmark###############################
### Benchmark between VSS, InterVar, Tapes (ref = ClinGen) ###
colnames(merge)

# Tapes_prob vss VSS2_score
benchmark_score <- merge %>%
  filter(
    Assertion %in% c(
      "Pathogenic",
      "Likely Pathogenic",
      "Benign",
      "Likely Benign"
    )
  ) %>%
  mutate(
    Truth = if_else(
      Assertion %in% c("Pathogenic", "Likely Pathogenic"),
      1, 0
    )
  ) %>%
  filter(
    !is.na(VSS2_score),
    !is.na(Tapes_prob)
  )

# Check sample size
nrow(benchmark_score)
table(benchmark_score$Truth)


## VSS ROC

roc_VSS <- roc(
  response = benchmark_score$Truth,
  predictor = benchmark_score$VSS2_score,
  levels = c(0, 1),
  direction = "<"
)

auc(roc_VSS)


## TAPES ROC

roc_Tapes <- roc(
  response = benchmark_score$Truth,
  predictor = benchmark_score$Tapes_prob,
  levels = c(0, 1),
  direction = "<"
)

auc(roc_Tapes)

benchmark_score %>%
  group_by(Truth) %>%
  summarise(
    n = n(),
    mean_Tapes_prob = mean(Tapes_prob),
    median_Tapes_prob = median(Tapes_prob)
  )

plot(
  roc_VSS,
  legacy.axes = TRUE,
  lwd = 2,
  main = "ROC Curve Comparison"
)

plot(
  roc_Tapes,
  add = TRUE,
  lwd = 2,
  lty = 2
)

abline(
  a = 0,
  b = 1,
  lty = 3
)

legend(
  "bottomright",
  legend = c(
    paste0("VSS (AUROC = ", round(auc(roc_VSS), 3), ")"),
    paste0("TAPES (AUROC = ", round(auc(roc_Tapes), 3), ")")
  ),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)

roc.test(
  roc_VSS,
  roc_Tapes,
  method = "delong",
  paired = TRUE
)

### ============================================================
### 3-class benchmark
### Truth: ClinGen
### Classes: P/LP, VUS, B/LB
### VSS rule:
###   score >= 6  -> P/LP
###   score 2-5   -> VUS
###   score <= 1  -> B/LB
### ============================================================

VSS_3class <- VSS %>%
  mutate(
    
    ## -------------------------
    ## ClinGen ground truth
    ## -------------------------
    Truth = case_when(
      
      Assertion %in% c(
        "Pathogenic",
        "Likely Pathogenic"
      ) ~ "P/LP",
      
      Assertion == "Uncertain Significance" ~ "VUS",
      
      Assertion %in% c(
        "Benign",
        "Likely Benign"
      ) ~ "B/LB",
      
      TRUE ~ NA_character_
    ),
    
    
    ## -------------------------
    ## VSS prediction
    ## -------------------------
    Prediction = case_when(
      
      VSS2_score >= 6 ~ "P/LP",
      
      VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
      
      VSS2_score <= 1 ~ "B/LB",
      
      TRUE ~ NA_character_
    )
    
  ) %>%
  
  ## remove unsupported ClinGen labels / missing score
  filter(
    !is.na(Truth),
    !is.na(Prediction)
  ) #12619 variants

table(
  Truth = VSS_3class$Truth,
  Prediction = VSS_3class$Prediction
)

class_levels <- c("P/LP", "VUS", "B/LB")

VSS_3class <- VSS_3class %>%
  mutate(
    Truth = factor(
      Truth,
      levels = class_levels
    ),
    
    Prediction = factor(
      Prediction,
      levels = class_levels
    )
  )

cm_3class <- confusionMatrix(
  data = VSS_3class$Prediction,
  reference = VSS_3class$Truth
)

cm_3class


### ============================================================
### 5-class benchmark
### Truth: ClinGen
### Classes: P, LP, VUS, LB, B
### ============================================================

VSS_5class <- VSS %>%
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
    Prediction = case_when(
      VSS2_score >= 10 ~ "P",
      VSS2_score >= 6 & VSS2_score <= 9 ~ "LP",
      VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
      VSS2_score >= 0 & VSS2_score <= 1 ~ "LB",
      VSS2_score < 0 ~ "B",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Truth),
    !is.na(Prediction)
  )


class_levels <- c("P", "LP", "VUS", "LB", "B")

VSS_5class <- VSS_5class %>%
  mutate(
    Truth = factor(Truth, levels = class_levels),
    Prediction = factor(Prediction, levels = class_levels)
  )

table(
  Prediction = VSS_5class$Prediction,
  Reference = VSS_5class$Truth
)

cm_5class <- confusionMatrix(
  data = VSS_5class$Prediction,
  reference = VSS_5class$Truth
)

cm_5class
