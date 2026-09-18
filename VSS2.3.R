### Required libraries
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(furrr)
library(purrr)
library(stringr)
library(reshape2)
library(optparse)
library(stringi)
library(data.table)
library(vcfR)

##### Setup cmd settings #####
###cmd: Rscript VSS2.2.R --data /path/to/your/snv.tables --ref /path/to/your/ref --output /path/to/results --BA1_restore

# Define the list of command-line options
option_list <- list(
  make_option(c("--data"), type = "character", help = "File path to the annotation files"),
  make_option(c("--ref"), type = "character", help = "File path to the reference tables"),
  make_option(c("--output"), type = "character", help = "File path to the output results"),
  make_option(c("--DP"), type = "double", default = 10, help = "Depth cutoff [default %default]", metavar = "Integer"),
  make_option(c("--AF"), type = "double", default = 0.05, help = "Rounded population allele frequency cutoff [default %default]", metavar = "Numerical"),
  make_option(c("--VSS2"), type = "double", default = 5, help = "VSS2 cutoff [default %default]", metavar = "Integer -2~13"),
  make_option(c("--QUAL_filter"), action = "store_true", default = FALSE, help = "QUAL filters (QD,SB,ReadPos) are on [default %default]"),
  make_option(c("--BA1_restore"), action = "store_true", default = FALSE, help = "BA1 variants are included [default %default]"),
  make_option(c("--phenotype"), type = "character", help = "Phenotype filter (comma-separated, no space)", metavar = "Phenotype keywords")
)

# Create OptionParser object
parser <- OptionParser(option_list = option_list)

# Parse the arguments
args <- parse_args(parser)

# Define required arguments and print error if they are missing
required_args <- c("data", "ref", "output")
missing_args <- required_args[!required_args %in% names(args) | is.null(args[required_args])]

if (length(missing_args) > 0) {
  print_help(parser)
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")), call. = FALSE)
}

# Assign required arguments
data_path <- args$data
ref_path <- args$ref
output_path <- args$output

# Assign optional arguments with defaults
DP_cutoff <- args$DP
AF_cutoff <- args$AF
VSS2_cutoff <- args$VSS2
QUAL_filter <- args$QUAL_filter  # TRUE if --QUAL_filter is provided
BA1_restore <- args$BA1_restore  # TRUE if --BA1_restore is provided
phenotype <- if (!is.null(args$phenotype)) {
  unlist(strsplit(args$phenotype, ","))
} else {
  NULL
}

# Print the parsed arguments
cat("----- Parsed Arguments -----\n")
cat("Data path:", data_path, "\n")
cat("Reference path:", ref_path, "\n")
cat("Output path:", output_path, "\n")

if (is.null(DP_cutoff)) {
  cat("DP cutoff:", 10, "\n")
} else {
  cat("DP cutoff:", DP_cutoff, "\n")
}

if (is.null(AF_cutoff)) {
  cat("Rounded AF cutoff:", 0.05, "\n")
} else {
  cat("Rounded AF cutoff:", AF_cutoff, "\n")
}

if (is.null(VSS2_cutoff)) {
  cat("VSS2 cutoff:", 5, "\n")
} else {
  cat("VSS2 cutoff:", VSS2_cutoff, "\n")
}

if (!is.null(phenotype)) {
  cat("Phenotype filter:", paste(phenotype, collapse = ", "), "\n")
} else {
  cat("Phenotype filter: Not specified\n")
}
cat("-----------------------------\n")


# Find VCF / VCF.GZ files in data_path
file_names <- list.files(
  path = data_path,
  pattern = "\\.vcf(\\.gz)?$",
  full.names = TRUE
)

cat("Input VCF:", basename(file_names))


#####----------Read reference files----------

csv_files <- list.files(
  ref_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

ref  <- NULL
wes  <- NULL
mask <- NULL


for (f in csv_files) {
  
  if (grepl("Used_columns_scores_vss", f)) {
    
    ref <- read.csv(f)
    cat("ref:", basename(f), "\n")
    
  } else if (grepl("Gene_panel_coordinates", f)) {
    
    wes <- read.csv(f)
    cat("panel:", basename(f), "\n")
    
  } else if (grepl("Common_variants_after_mask", f)) {
    
    mask <- read.csv(f)
    cat("mask:", basename(f), "\n")
  }
}


#####----------Save R.workspace----------
save_workspace <- function() {
  
  tryCatch({
    
    save.image(
      file = file.path(output_path, "workspace.RData")
    )
    
    cat(
      "Workspace saved to:",
      file.path(output_path, "workspace.RData"),
      "\n"
    )
    
  }, error = function(e) {
    
    cat("Error while saving workspace\n")
    print(e)
    
  })
}

#####----------Read vcfs by vcfR----------
process_vcf <- function(vcf_file){
  
  cat("Processing:",
      basename(vcf_file),
      "\n")
  
  # read vcf
  vcf <- tryCatch({
    
    read.vcfR(vcf_file)
    
  }, error=function(e){
    
    cat("Cannot read:",
        basename(vcf_file),
        "\n")
    
    print(e)
    
    return(NULL)
  })
  
  if(is.null(vcf)){
    return(NULL)
  }
  
  # Change to tidy dataframe
  tidy_all <- vcfR2tidy(
    vcf,
    info_only=FALSE)
  
  gt_df <- tidy_all$gt
  fix_df  <- tidy_all$fix
  
  data <- cbind(gt_df,fix_df)
  data <- data[, !duplicated(colnames(data))]
  
  return(data)
}


table_list <- lapply(file_names, process_vcf)
table_list <- table_list[!sapply(table_list, is.null)]

names(table_list) <- sapply(table_list, function(df) {
  ind <- unique(df$Indiv)
  if (length(ind) == 1) ind else paste(ind, collapse = "_")
})

table_list <- lapply(table_list, function(df) {
  if (is.null(df)) return(NULL)
  
  # 移除 ChromKey, FILTER, DB
  df <- df[, !colnames(df) %in% c("ChromKey", "FILTER", "DB", "InbreedingCoeff"), drop = FALSE]
  
  # 欄位重排
  df <- df[, c(
    "CHROM", "POS", "REF", "ALT", "ID",
    setdiff(colnames(df), c("CHROM", "POS", "REF", "ALT", "ID"))
  )]
  
  return(df)
})



# Examine duplicate rows in the data (If there is warning messages means there are error in the input file)
check_duplicates <- function(df, name){
  
  dup_variant <- df %>%
    count(CHROM, POS, REF, ALT) %>%
    filter(n > 1)
  
  if(nrow(dup_variant) > 0){
    
    warning(
      paste0(
        nrow(dup_variant),
        " duplicated variants found in ",
        name
      )
    )
    
    print(dup_variant)
    
  } else {
    
    cat(
      "No duplicated variants:",
      name,
      "\n"
    )
  }
}

for(i in seq_along(table_list)){
  
  check_duplicates(
    table_list[[i]],
    basename(file_names[i])
  )
}


# Wrap your main script logic in a tryCatch block
tryCatch({
  
  ##### Start of VSS2.3 main script  ##### 
  
  ##### ---------- Data reformatting---------- 
  
  ## 清理x3b(;)跟x3d(=)跟., 並針對dbNSFP_OMIM_id欄位做6位數逗號格式化
  clean_raw_data <- function(data) {
    data %>%
      # Step 1：factor → character
      mutate(across(where(is.factor), as.character)) %>%
      
      # Step 2：清理特殊字元
      mutate(across(where(is.character), ~ {
        
        s <- stri_replace_all_fixed(., "\\x3b", ";")  # \x3b → ";"
        s <- stri_replace_all_fixed(s, "\\x3d", "=")  # \x3d → "="
        s <- stri_replace_all_fixed(s, ".,", "")       # 移除 ".,"
        
        # split by ";"
        parts <- stri_split_fixed(s, ";", simplify = FALSE)
        
        cleaned <- lapply(parts, function(x) {
          x <- stri_trim_both(x)
          x <- x[!(is.na(x) | x %in% c(".", ""))]
          
          if (length(x) == 0) return(NA_character_)
          paste(x, collapse = ";")
        })
        
        unlist(cleaned)
      })) %>%
      
      # Step 3：OMIM formatting變6位數字
      mutate(
        dbNSFP_OMIM_id = {
          tmp <- dbNSFP_OMIM_id
          
          non_na_idx <- !is.na(tmp)
          
          tmp[non_na_idx] <- sapply(
            strsplit(tmp[non_na_idx], ","),
            function(x) {
              paste(sprintf("%06d", as.numeric(trimws(x))), collapse = ";")
            }
          )
          
          tmp
        }
      )
  }
  
  
  ## 清理重複值
  dedup_semicolon_cells <- function(data) {
    data %>%
      mutate(across(where(is.character), ~ {
        s <- .                              
        # 1) 用 stringi 在 C 層級把每個 cell 用 ";" 拆成 list
        parts_list <- stri_split_fixed(s, ";", omit_empty = FALSE)
        
        # 2) 對每個拆好的元素做去 trim、移除空字串、去重複、再 paste 回字串
        out <- vapply(parts_list, FUN.VALUE = character(1), USE.NAMES = FALSE, FUN = function(parts) {
          # stringi 的 split 在遇到 NA 原始值會回傳 NA（非向量），所以先處理 NA
          if (length(parts) == 1 && is.na(parts)) return(NA_character_)
          
          # trim 兩端空白（在 C 層面執行）
          parts <- stri_trim_both(parts)
          # 移除完全是空字串
          parts <- parts[parts != ""]
          # 移除完全等於 NA（避免 weird cases）
          parts <- parts[!is.na(parts)]
          # 若完全沒有剩下，就回 NA
          if (length(parts) == 0) return(NA_character_)
          # 保留第一個出現順序的 unique
          parts <- unique(parts)
          # 用 stringi join 回字串（C 層級）
          stri_join(parts, collapse = ";")
        })
        
        out
      }))
  }
  
  ## Merge columns
  merge_columns <- function(data) {
    data$phastCons100way_vertebrate <- ifelse(is.na(data$phastCons100way_vertebrate), data$verPhCons, data$phastCons100way_vertebrate)
    data$phyloP100way_vertebrate <- ifelse(is.na(data$phyloP100way_vertebrate), data$verPhyloP, data$phyloP100way_vertebrate)
    data$CADD_phred <- ifelse(is.na(data$CADD_phred), data$PHRED, data$CADD_phred)
    data <- data %>% select(-verPhCons, -verPhyloP, -PHRED)
    
    return(data)
  }
  
  ## Add gt_AF
  add_gt_AF <- function(df, ad_col = "gt_AD") {
    
    df$gt_AF <- vapply(
      strsplit(as.character(df[[ad_col]]), ","),
      function(x) {
        if(length(x) >= 2) {
          ref <- as.numeric(x[1])
          alt <- as.numeric(x[2])
          
          alt / (ref + alt)
        } else {
          NA_real_
        }
      },
      numeric(1)
    )
    
    return(df)
  }
  
  ## Assign WES/ECS/INF gene panels
  check_panel_hits <- function(data, wes) {
    
    # 建立結果欄位，預設 0
    data$WES <- 0
    data$ECS <- 0
    data$INF <- 0
    data$ECS680 <- 0
    data$Paralogous <- 0
    
    # 逐筆 variant 判斷
    for (i in seq_len(nrow(data))) {
      chr <- data$CHROM[i]
      pos <- data$POS[i]
      
      # 抓出同染色體的區間
      wes_chr <- wes[wes$CHROM == chr, ]
      
      if (nrow(wes_chr) > 0) {
        # 找出落在 Start ~ End 的 row（可能有多筆）
        hit_rows <- wes_chr$Start <= pos & wes_chr$End >= pos
        
        if (any(hit_rows)) {
          # 如果落在區間內，把 panel hit 填上
          data$WES[i] <- as.numeric(any(wes_chr$WES[hit_rows] == 1))
          data$ECS[i] <- as.numeric(any(wes_chr$ECS[hit_rows] == 1))
          data$INF[i] <- as.numeric(any(wes_chr$INF[hit_rows] == 1))
          data$ECS680[i] <- as.numeric(any(wes_chr$ECS680[hit_rows] == 1))
          data$Paralogous[i] <- as.numeric(any(wes_chr$Paralogous[hit_rows] == 1))
        }
      }
    }
    
    return(data)
  }
  
  
  # Apply the rename and replace functions to each data frame in the list
  plan(multisession)
  
  table_list <- future_map(
    table_list,
    function(df) {
      df %<>% clean_raw_data()
      df %<>% dedup_semicolon_cells()
      df %<>% add_gt_AF()
      df %<>% merge_columns()
      
      df %<>% mutate(across(
        where(is.character),
        ~ na_if(., "") |> na_if("NA") |> na_if("N/A") |> na_if("NULL")
      ))
      
      df %<>% check_panel_hits(wes = wes)   # ← 正確傳參數
      
      df %<>% mutate(
        Mask_common = as.integer(
          paste(CHROM, POS, REF, ALT, sep = "_") %in%
            paste(mask$CHROM, mask$POS, mask$REF, mask$ALT, sep = "_")
        )
      )
      
      return(df)
    }
  )
  
  
  ##### ----------Step 1. Select any pop AF <= 0.05 & BA1-excluded variants by CLNID ----------
  
  # 指定要看的AF欄位
  af_cols <- c(
    "TWB_AF_1517",
    "gnomad_exomes_AF",
    "gnomad_exomes_AF_eas", 
    "gnomAD4.1_joint_AF",
    "gnomAD4.1_joint_EAS_AF",
    "AllofUs_ALL_AF", 
    "AllofUs_EAS_AF",
    "RegeneronME_ALL_AF",
    "RegeneronME_EAS_AF",
    "ALFA_Total_AF",
    "ALFA_East_Asian_AF") #共11欄位
  
  PM2_cols <- c("TWB_AF_1517",
                "gnomad_exomes_AF",
                "gnomad_exomes_AF_eas")
  
  
  # BA1 excluded pathogenic variants
  exclude_CLNID <- c("1018", "17023", "10", "9", "2551", "2552", "217689", "3830", "1900")
  
  ACMG_AF <- function(data) {
    
    # Convert relevant columns to numeric
    data <- data %>%
      mutate(across(all_of(af_cols), ~ suppressWarnings(as.numeric(.))))
    
    # Add ACMG_AF (PM2 and BA1)
    data <- data %>%
      mutate(ACMG_AF = case_when(
        
        # PM2: 只看gnomad_exomes_AF | gnomad_exomes_AF_eas <= 1e-4 或 加上TWB全為NA者
        (if_all(all_of(PM2_cols), is.na)) |
          (gnomad_exomes_AF <= 0.0001) |
          (gnomad_exomes_AF_eas <= 0.0001) ~ "PM2",
        
        # BA1: CLNID 不在排除清單，且任一 AF > 0.05
        !CLNID %in% exclude_CLNID & if_any(all_of(PM2_cols), ~ . > 0.05) ~ "BA1",
        
        # 其他情況
        TRUE ~ NA_character_
      ))
    
    return(data)
  }
  
  plan(multisession)
  data1 <- future_map(table_list, ACMG_AF)
  
  
  # Generate a report of the AF_filtered variant no.
  data1_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF)))
  )
  
  print(data1_rows)
  
  
  ##### ----------Step 2. Mutation type + Clinvar evidence----------
  
  ACMG_MutType <- function(data) {
    
    # --- Step 1: 初步 assign ACMG_MutType ---
    data2 <- data %>%
      mutate(ACMG_MutType = case_when(
        
        # PVS1
        (ExonicFunc.ensGene %in% c(
          "frameshift_insertion", "frameshift_deletion", "frameshift_substitution",
          "frameshift_block_substitution", "stopgain", "startloss", "startgain"
        ) |
          ExonicFunc.refGene %in% c(
            "frameshift_insertion", "frameshift_deletion", "frameshift_substitution",
            "frameshift_block_substitution", "stopgain", "startloss", "startgain"
          )) ~ "PVS1",
        
        # PM4
        (ExonicFunc.ensGene %in% c("stoploss", "nonframeshift_insertion", "nonframeshift_deletion") |
           ExonicFunc.refGene %in% c("stoploss", "nonframeshift_insertion", "nonframeshift_deletion")) ~ "PM4",
        
        # PS1: missense + ClinVar pathogenic (no conflicting)
        ((ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) &
           grepl("athogenic", CLNSIG) & !grepl("Conflicting", CLNSIG)) |
          ((ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) &
             grepl("athogenic", CLNSIG) & !grepl("Conflicting", CLNSIG)) ~ "PS1",
        
        # PM5: missense + ClinVar VUS
        ((ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) &
           grepl("Uncertain", CLNSIG)) |
          ((ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) &
             grepl("Uncertain", CLNSIG)) ~ "PM5",
        
        # BP7: Synonymous
        (ExonicFunc.ensGene == "synonymous_SNV" | ExonicFunc.refGene == "synonymous_SNV") & !grepl("athogenic", CLNSIG) ~ "BP7",
        
        # BP6: Benign missense
        (((ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution",
                                     "nonframeshift_block_substitution", "unknown")) | is.na(ExonicFunc.ensGene)) &
           grepl("enign", CLNSIG)) |
          (((ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution",
                                       "nonframeshift_block_substitution", "unknown")) | is.na(ExonicFunc.refGene)) &
             grepl("enign", CLNSIG)) ~ "BP6",
        
        # PP2: missense + ClinVar no data OR CLNSIGCONF has no pathogenic evidence
        ((ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) |
           (ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution"))) &
          (is.na(CLNSIG) | (grepl("Conflicting", CLNSIG) & !grepl("athogenic", CLNSIGCONF))) ~ "PP2",
        
        # PP5: no mutation type but ClinVar pathogenic (no conflicting)
        (is.na(ExonicFunc.ensGene) & is.na(ExonicFunc.refGene) & grepl("athogenic", CLNSIG) & !grepl("Conflicting", CLNSIG)) ~ "PP5",
        
        TRUE ~ NA_character_
      ))
    
    # --- Step 2: 定義函數處理 CLNSIGCONF ---
    calculate_evidence_counts <- function(clnsigconf) {
      items <- unlist(strsplit(clnsigconf, "\\|"))
      
      pathogenic_count <- 0
      benign_count <- 0
      
      for (item in items) {
        if (grepl("athogenic", item)) {
          count <- as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", item))
          pathogenic_count <- pathogenic_count + count
        } else if (grepl("enign", item)) {
          count <- as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", item))
          benign_count <- benign_count + count
        }
      }
      
      list(pathogenic_count = pathogenic_count, benign_count = benign_count)
    }
    
    # --- Step 3: 再用 CLNSIGCONF 判斷 ---
    data_filtered <- data2 %>%
      filter(
        ((is.na(ExonicFunc.ensGene) | is.na(ExonicFunc.refGene)) & !is.na(CLNSIGCONF)) |
          (ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution") & !is.na(CLNSIGCONF)) |
          (ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution") & !is.na(CLNSIGCONF))
      ) %>%
      rowwise() %>%
      mutate(
        ACMG_MutType = {
          conf_list <- calculate_evidence_counts(CLNSIGCONF)
          
          if (
            ((ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution")) |
             (ExonicFunc.refGene %in% c("nonsynonymous_SNV", "nonframeshift_substitution", "nonframeshift_block_substitution"))) &
            grepl("Conflicting", CLNSIG)
          ) {
            if (conf_list$pathogenic_count > conf_list$benign_count) {
              "PS1"
            } else {
              "PM5"
            }
          } else if (
            (is.na(ExonicFunc.ensGene) & is.na(ExonicFunc.refGene)) & !is.na(CLNSIGCONF)
          ) {
            if (conf_list$pathogenic_count > conf_list$benign_count) {
              "PP5"
            } else {
              NA_character_
            }
          } else {
            NA_character_
          }
        }
      ) 
    
    # --- Step 4: 合併結果 ---
    data2 <- data2 %>%
      left_join(
        data_filtered %>% 
          select(CHROM, POS, REF, ALT, ACMG_MutType) %>%
          filter(!is.na(ACMG_MutType)), 
        by = c("CHROM", "POS", "REF", "ALT")
      ) %>%
      mutate(ACMG_MutType = coalesce(ACMG_MutType.x, ACMG_MutType.y)) %>%
      select(-ACMG_MutType.x, -ACMG_MutType.y) 
    
    return(data2)
  }
  
  
  plan(multisession)
  data2 <- future_map(data1, ACMG_MutType)
  
  
  # Generate a report of the ACMG_MutType
  data2_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF))),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType)))
  )
  
  print(data2_rows)
  
  
  #####----------Step 3. Splicing evidence----------
  ACMG_Splicing <- function(data, ref, sp_coltype = "Splicing_predictor") {
    
    # 取出 splicing predictor 欄位
    sp_cols <- ref %>% 
      filter(Coltype == sp_coltype) %>% 
      pull(Colname)
    
    # 定義單一值對照 cutoff 的 helper
    assign_splice_level <- function(value, colname) {
      if (is.na(value)) return(NA_character_)
      
      cutoffs <- ref %>%
        filter(Colname == colname) %>%
        select(Supporting_pathogenic, Moderate_pathogenic, Strong_pathogenic) %>%
        unlist()
      
      if (!is.na(cutoffs["Strong_pathogenic"]) && value >= cutoffs["Strong_pathogenic"]) {
        return("PS_S")
      } else if (!is.na(cutoffs["Moderate_pathogenic"]) && value >= cutoffs["Moderate_pathogenic"]) {
        return("PM_S")
      } else if (!is.na(cutoffs["Supporting_pathogenic"]) && value >= cutoffs["Supporting_pathogenic"]) {
        return("PP_S")
      } else {
        return(NA_character_)
      }
    }
    
    # 先計算每個 predictor 的等級矩陣
    splice_matrix <- sapply(sp_cols, function(col) {
      v <- data[[col]]
      sapply(seq_along(v), function(i) assign_splice_level(v[i], col))
    })
    
    # 定義等級 ranking
    level_rank <- c("PP_S" = 1, "PM_S" = 2, "PS_S" = 3)
    
    # 取每行最高等級
    acmg <- apply(splice_matrix, 1, function(lv) {
      lv <- lv[!is.na(lv)]
      if(length(lv) == 0) return(NA_character_)
      names(which.max(level_rank[lv]))
    })
    
    # 加上 canonical splice 判斷，直接覆蓋為 PVS1
    acmg <- ifelse(
      data$Func.ensGene == "splicing" | data$Func.refGene == "splicing" |
        data$Consequence %in% c("CANONICAL_SPLICE", "SPLICE_SITE"),
      "PVS1",
      acmg
    )
    
    # 回傳原表加上 ACMG_splicing 欄位
    data$ACMG_Splicing <- acmg
    return(data)
  }
  
  plan(multisession)
  data3 <- future_map(data2, ~ACMG_Splicing(.x, ref, sp_coltype = "Splicing_predictor"))
  
  
  # Generate a report of the ACMG_Splicing
  data3_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF))),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType))),
    ACMG_Splicing = sapply(data3, function(df) sum(!is.na(df$ACMG_Splicing)))
  )
  
  print(data3_rows)
  
  
  #####----------Step 4. Inter-protein domain evidence ----------
  
  ACMG_InterProDomain <- function(data3){
    
    data4 <- data3 %>%
      mutate(ACMG_InterProDomain = ifelse(!is.na(Interpro_domain) | !is.na(VEP_domain) , "PM1", NA))
    
    return(data4)  
  }  
  
  plan(multisession)
  data4 <- future_map(data3, ACMG_InterProDomain)
  
  # Generate a report of the ACMG_InterProDomain
  data4_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF))),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType))),
    ACMG_Splicing = sapply(data3, function(df) sum(!is.na(df$ACMG_Splicing))),
    ACMG_InterProDomain = sapply(data4, function(df) sum(!is.na(df$ACMG_InterProDomain)))
  )
  
  print(data4_rows)
  
  
  #####----------Step 5. In-silico predictor evidence----------
  cols <- ref %>% filter(Coltype %in% c("Insilico_predictor", "Conservation_predictor")) %>% pull(Colname)
  
  # Convert to numeric with split ; and take max
  convert_to_numeric <- function(df, cols) {
    df %>%
      select(all_of(cols)) %>%
      mutate(across(
        all_of(cols),
        ~ suppressWarnings(
          sapply(strsplit(as.character(.x), ";"), function(v) {
            v_num <- suppressWarnings(as.numeric(v))
            if (all(is.na(v_num))) {
              NA_real_
            } else {
              v_num[which.max(abs(v_num))]   # 取絕對值最大的原始數字
            }
          })
        )
      )) %>%
      bind_cols(select(df, CHROM, POS, REF, ALT))
  }
  
  data5 <- map(data4, ~ convert_to_numeric(.x, cols))
  
  
  # Five scores have to be inverted
  data5 <- map(data5, ~ {
    .x %>%
      mutate(
        ESM1b_score = -1 * ESM1b_score
      )
  })
  
  ## Assign score range and interpretation to each tool
  assign_category <- function(value, cutoffs) {
    if (is.na(value)) {
      return(NA)
    }
    if (!is.na(cutoffs$Strong_benign) && value <= cutoffs$Strong_benign) {
      return("Strong_benign")
    } else if (!is.na(cutoffs$Moderate_benign) && value <= cutoffs$Moderate_benign) {
      return("Moderate_benign")
    } else if (!is.na(cutoffs$Supporting_benign) && value <= cutoffs$Supporting_benign) {
      return("Supporting_benign")
    } else if (!is.na(cutoffs$Strong_pathogenic) && value >= cutoffs$Strong_pathogenic) {
      return("Strong_pathogenic")
    } else if (!is.na(cutoffs$Moderate_pathogenic) && value >= cutoffs$Moderate_pathogenic) {
      return("Moderate_pathogenic")
    } else if (!is.na(cutoffs$Supporting_pathogenic) && value >= cutoffs$Supporting_pathogenic) {
      return("Supporting_pathogenic")
    } else {
      return("VUS")
    }
  }
  
  score <- ref %>% filter(Coltype %in% c("Insilico_predictor", "Conservation_predictor"))
  assign_categories <- function(df, score) {
    for (col in cols) {
      if (col %in% names(df)) {
        cutoffs <- score %>% filter(Colname == col)
        if (nrow(cutoffs) == 1) {
          df <- df %>%
            mutate("{col}_category" := sapply(get(col), assign_category, cutoffs = cutoffs))
        }
      }
    }
    df
  }
  
  # Apply category assignment to each dataframe in data5
  data5 <- map(data5, ~ assign_categories(.x, score))
  
  # Assign category weights
  weights <- c(
    "Strong_benign" = -4,
    "Moderate_benign" = -2,
    "Supporting_benign" = -1,
    "VUS" = 0,
    "Supporting_pathogenic" = 1,
    "Moderate_pathogenic" = 2,
    "Strong_pathogenic" = 4
  )
  
  # Calculate_weighted_mean
  calculate_weighted_mean <- function(df_row, weights) {
    # Select colnames with "_category"
    df_row <- df_row[grep("_category", names(df_row))]
    
    # Remove NA and sum weight*count of category
    counts <- sapply(df_row, function(column) {
      if (all(is.na(column))) {
        return(NA)
      }
      tab <- table(column)
      if (length(tab) == 0) {
        return(NA)
      }
      sum(tab * weights[names(tab)])
    })
    
    # Calculate weighted mean
    if (all(is.na(counts))) {
      weighted_mean <- NA
    } else {
      weighted_mean <- sum(counts, na.rm = TRUE) / sum(!is.na(df_row))
    }
    
    return(weighted_mean)
  }
  
  
  # Apply calculate_weighted_mean function to each dataframe in table_list
  data5 <- map(data5, ~ {
    df <- .x
    df$Weighted_mean <- apply(df, 1, calculate_weighted_mean, weights = weights)
    # Assign ACMG_predictors based on Weighted_mean
    df$ACMG_Predictors <- ifelse(df$Weighted_mean >= 0, "PP3", "BP4")
    df
  })
  
  
  
  # Generate a report of the ACMG_Predictors
  data5_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(df$ACMG_AF == "PM2", na.rm = TRUE)),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType))),
    ACMG_Splicing = sapply(data3, function(df) sum(!is.na(df$ACMG_Splicing))),
    ACMG_InterProDomain = sapply(data4, function(df) sum(!is.na(df$ACMG_InterProDomain))),
    ACMG_Predictors = sapply(data5, function(df) sum(!is.na(df$ACMG_Predictors)))
  )
  
  print(data5_rows)
  
  
  # Select data5's Weighted_mean & ACMG_predictors，merge back to data4 by CHROM POS REF ALT
  data5_processed <- lapply(data5, function(df) {
    df %>% 
      select(CHROM, POS, REF, ALT, Weighted_mean, ACMG_Predictors)
  })
  
  final_data_list <- map2(data4, data5_processed, ~ inner_join(.x, .y, by = c("CHROM", "POS", "REF", "ALT")))
  
  
  #####----------Step6. Calculate final VSS2 score and make interpretation (BA1 not counted)---------- 
  
  calculate_VSS2 <- function(df) {
    # 權重定義
    weights <- c("PVS" = 8, "PS" = 4, "PM" = 2, "PP" = 1, "BP" = -1)
    
    # 選 ACMG 欄位
    acmg_cols <- c("ACMG_AF", "ACMG_MutType", "ACMG_Splicing", "ACMG_InterProDomain", "ACMG_Predictors")
    acmg_columns <- df %>% select(all_of(acmg_cols))
    
    # 幫助函數：依部分匹配給權重
    assign_weight <- function(value) {
      # 對 NA 或空字串回傳 NA（讓我們可以分辨「真正沒有證據」與「有證據且是0/負分」）
      if (is.na(value) || trimws(value) == "") return(NA_real_)
      for (k in names(weights)) {
        if (grepl(k, value, ignore.case = TRUE)) return(weights[[k]])
      }
      # 如果字串裡面沒有找到任何 PVS/PS/PM/PP/BP，就回傳 0（代表沒有那類證據）
      return(0)
    }
    
    # 計算每個欄位權重（保留 NA）
    acmg_weighted <- acmg_columns %>%
      mutate(across(everything(), ~ sapply(., assign_weight), .names = "weighted_{col}"))
    
    # MutType 與 Splicing 的合併規則：
    # - 若只有一個非 NA，就採那個（包含負分）
    # - 若兩個都有，取較大值
    # - 若兩個都 NA，回傳 0（以免影響後面加總）
    acmg_weighted <- acmg_weighted %>%
      mutate(
        weighted_Mut_Splice = case_when(
          !is.na(weighted_ACMG_MutType) & is.na(weighted_ACMG_Splicing) ~ weighted_ACMG_MutType,
          is.na(weighted_ACMG_MutType) & !is.na(weighted_ACMG_Splicing) ~ weighted_ACMG_Splicing,
          !is.na(weighted_ACMG_MutType) & !is.na(weighted_ACMG_Splicing) ~ pmax(weighted_ACMG_MutType, weighted_ACMG_Splicing),
          TRUE ~ 0
        ),
        Mut_Splice_source = case_when(
          !is.na(weighted_ACMG_MutType) & is.na(weighted_ACMG_Splicing) ~ "MutType",
          is.na(weighted_ACMG_MutType) & !is.na(weighted_ACMG_Splicing) ~ "Splicing",
          !is.na(weighted_ACMG_MutType) & !is.na(weighted_ACMG_Splicing) & weighted_ACMG_Splicing >= weighted_ACMG_MutType ~ "Splicing",
          !is.na(weighted_ACMG_MutType) & !is.na(weighted_ACMG_Splicing) & weighted_ACMG_MutType > weighted_ACMG_Splicing ~ "MutType",
          TRUE ~ NA_character_
        )
      )
    
    # 計算 PVS/PS/PM/PP/BP count（維持您原本的邏輯）
    count_acmg <- acmg_columns %>%
      mutate(
        Mut_Splice_value = mapply(function(mut, splice, source){
          if (is.na(source)) return(NA_character_)
          if (source == "Splicing") return(splice)
          else return(mut)
        }, ACMG_MutType, ACMG_Splicing, acmg_weighted$Mut_Splice_source),
        
        PVS_count = ifelse(grepl("PVS", Mut_Splice_value, ignore.case = TRUE), 1, 0),
        PS_count  = ifelse(grepl("PS", Mut_Splice_value, ignore.case = TRUE), 1, 0),
        PM_count  = ifelse(grepl("PM", Mut_Splice_value, ignore.case = TRUE), 1, 0) +
          ifelse(grepl("PM", ACMG_AF, ignore.case = TRUE), 1, 0) +
          ifelse(grepl("PM", ACMG_InterProDomain, ignore.case = TRUE), 1, 0),
        PP_count  = ifelse(grepl("PP", Mut_Splice_value, ignore.case = TRUE), 1, 0) +
          ifelse(grepl("PP", ACMG_Predictors, ignore.case = TRUE), 1, 0),
        BP_count  = ifelse(grepl("BP", Mut_Splice_value, ignore.case = TRUE), 1, 0) +
          ifelse(grepl("BP", ACMG_Predictors, ignore.case = TRUE), 1, 0)
      )
    
    # 計算 VSS2（維持您原有的 NA 判斷）
    other_weights_sum <- rowSums(acmg_weighted %>%
                                   select(starts_with("weighted_ACMG_")) %>%
                                   select(-weighted_ACMG_MutType, -weighted_ACMG_Splicing), # 這裡只加 AF / InterProDomain / Predictors
                                 na.rm = TRUE)
    
    # 注意：為了確保上面選的欄位正確，明確列出要加總的 weighted 欄位也可以：
    # other_weights_sum <- rowSums(acmg_weighted %>% select(weighted_ACMG_AF, weighted_ACMG_InterProDomain, weighted_ACMG_Predictors), na.rm = TRUE)
    
    all_acmg_na <- rowSums(is.na(acmg_columns)) == ncol(acmg_columns)
    all_other_acmg_na <- rowSums(is.na(acmg_columns %>% select(-ACMG_AF))) == (ncol(acmg_columns) - 1)
    acmg_af_ba1 <- acmg_columns$ACMG_AF == "BA1"
    
    df <- df %>%
      mutate(
        VSS2_score = ifelse(all_acmg_na | (acmg_af_ba1 & all_other_acmg_na), NA,
                            other_weights_sum + acmg_weighted$weighted_Mut_Splice),
        PVS_count = count_acmg$PVS_count,
        PS_count  = count_acmg$PS_count,
        PM_count  = count_acmg$PM_count,
        PP_count  = count_acmg$PP_count,
        BP_count  = count_acmg$BP_count
      )
    
    return(df)
  }
  
  
  # Apply the function to each dataframe in the list
  plan(multisession)
  data_final_all <- future_map(final_data_list, calculate_VSS2)
  
  
  #####---------Step 7. Customized filtering--------------
  
  
  ## Function to call Inheritance mode
  extract_Inheritance <- function(df) {
    
    # 定義 inheritance 對照表
    Inheritance_map <- c(
      "Autosomal_recessive_inheritance"   = "AR",
      "Autosomal_dominant_inheritance"    = "AD",
      "X-linked_inheritance"              = "XL",
      "X-linked_dominant_inheritance"     = "XLD",
      "X-linked_recessive_inheritance"    = "XLR",
      "Polygenic_inheritance"             = "PG",
      "Digenic_inheritance"               = "DG",
      "Oligogenic_inheritance"            = "OG",
      "Mitochondrial_inheritance"         = "Mito",
      "Non-Mendelian_inheritance"         = "Non-Men"
    )
    
    df <- df %>%
      rowwise() %>%
      mutate(
        # 萃取 inheritance terms
        Inheritance_term = {
          terms <- unlist(str_split(dbNSFP_HPO_name, "[,;]"))  # , 或 ; 都拆開terms 
          # 去掉開頭的多個 . 和前後空格
          terms <- str_replace_all(terms, "^[\\.,\\s]+", "")  # 去掉開頭的 . 或 , 或空白
          terms <- str_replace_all(terms, "^\\s+|\\s+$", "")  # 去掉前後空白
          # 篩選出 _inheritance 結尾
          inh <- terms[str_detect(terms, regex("_inheritance$", ignore_case = TRUE))]
          inh <- inh[!inh %in% c("", ".", NA)]   # 過濾空值或 "."
          if (length(inh) == 0) NA_character_ else paste(unique(inh), collapse = ", ")
        },
        
        # 轉換成縮寫
        Inheritance = {
          if (is.na(Inheritance_term)) {
            NA_character_
          } else {
            replaced <- str_replace_all(Inheritance_term, Inheritance_map)
            abbr <- str_split(replaced, ", ")[[1]] %>% .[!. %in% c("", ".", NA)]
            if (any(str_detect(abbr, "_inheritance"))) NA_character_ else paste(sort(unique(abbr)), collapse = ", ")
          }
        }
      ) %>%
      ungroup() %>%
      # 將 NA 統一成 "NA" 字串
      mutate(
        Inheritance_term = ifelse(is.na(Inheritance_term), "NA", Inheritance_term),
        Inheritance = ifelse(is.na(Inheritance), "NA", Inheritance)
      )
    
    return(df)
  }
  
  
  plan(multisession)
  data_final_all <- future_map(data_final_all, ~extract_Inheritance(.x))
  
  
  ## Function to map phenotype keywords 
  map_phenotype <- function(df, phenotype = NULL) {
    
    if (is.null(phenotype) || length(phenotype) == 0) return(df)
    
    phenotype <- unlist(str_split(phenotype, ","))
    phenotype <- str_trim(phenotype)
    phenotype <- phenotype[!is.na(phenotype) & phenotype != ""]
    
    hpo_list <- str_split(df$dbNSFP_HPO_name, ";")
    
    phenotype_match <- lapply(hpo_list, function(x) {
      
      x <- str_trim(x)
      x <- x[
        !is.na(x) &
          x != "" &
          toupper(x) != "NA"
      ]
      
      if (length(x) == 0) {
        return(character(0))
      }
      
      hit <- sapply(x, function(term) {
        any(
          sapply(
            phenotype,
            function(p) {
              isTRUE(
                str_detect(
                  term,
                  fixed(p, ignore_case = TRUE)
                )
              )
            }
          )
        )
      })
      
      # 只保留 unique matched phenotype terms
      unique(x[hit])
    })
    
    df %>%
      mutate(
        Phenotype_match = sapply(
          phenotype_match,
          function(x) {
            if (length(x) == 0) {
              NA_character_
            } else {
              paste(x, collapse = ";")
            }
          }
        ),
        Phenotype_hit = lengths(phenotype_match)
      ) %>%
      filter(Phenotype_hit > 0)
  }
  
  
  ## Adding QD_filter, SB_filter, ReadPos_filter
  data_final_all <- future_map(
    data_final_all,
    ~ .x %>%
      mutate(
        QD_filter = if_else(QD < 2, "Fail", "Pass"),
        SB_filter = if_else(FS > 60 | SOR > 3, "Fail", "Pass"),
        ReadPos_filter = if_else(
          is.na(ReadPosRankSum) | ReadPosRankSum >= -8,
          "Pass",
          "Fail"
        )
      )
  )
  
  
  ## Customized filter
  customized_filter <- function(df, 
                                DP_cutoff = 10, 
                                AF_cutoff = 0.05, 
                                VSS2_cutoff = 5, 
                                BA1 = TRUE, 
                                phenotype = NULL,
                                QUAL_filter = FALSE) {
    
    # Apply DP cutoff
    df <- df %>% 
      filter(gt_DP >= DP_cutoff)
    
    
    # Apply AF cutoff
    df <- df %>%
      filter(
        CLNID %in% exclude_CLNID | 
          if_all(all_of(af_cols), is.na) |
          if_any(all_of(af_cols), ~ round(., 2) <= AF_cutoff)
      )
    
    # Apply VSS2 cutoff
    df <- df %>% 
      filter(VSS2_score >= VSS2_cutoff)
    
    # Exclude BA1 variants unless BA1 is TRUE
    if (!BA1) {
      df <- df %>% 
        filter(is.na(ACMG_AF) | ACMG_AF != "BA1")
    }
    
    # Apply QUAL filters
    if (QUAL_filter) {
      df <- df %>%
        filter(
          QD_filter == "Pass",
          SB_filter == "Pass",
          ReadPos_filter == "Pass" | is.na(ReadPos_filter)
        )
    }
    
    # Apply phenotype filtering
    if (!is.null(phenotype)) {
      df <- map_phenotype(df, phenotype)
    }
    
    return(df)
  }
  
  # Apply the function to each dataframe in the list
  plan(multisession)
  data_final_filtered <- future_map(data_final_all, 
                                    ~customized_filter(.x, DP = DP_cutoff, AF_cutoff = AF_cutoff, VSS2_cutoff = VSS2_cutoff, BA1 = BA1_restore, 
                                                       QUAL_filter = QUAL_filter, phenotype = phenotype))
  
  
  ## Generate a report for the filtered variants
  data_final_filtered_rows <- data.frame(
    Original = sapply(table_list, nrow),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF))),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType))),
    ACMG_Splicing = sapply(data3, function(df) sum(!is.na(df$ACMG_Splicing))),
    ACMG_InterProDomain = sapply(data4, function(df) sum(!is.na(df$ACMG_InterProDomain))),
    ACMG_Predictors = sapply(data5, function(df) sum(!is.na(df$ACMG_Predictors))),
    Variants_Filtered = sapply(data_final_filtered, function(df) nrow(df)))
  
  print(data_final_filtered_rows)
  
  
  #####----------Step 8. Output result files----------
  
  # Add ACMG combination & VSS2_classification/version
  data_final_filtered <- future_map(data_final_filtered, ~ .x %>%
                                      mutate(
                                        ACMG = paste(ACMG_AF, ACMG_MutType, ACMG_Splicing, ACMG_InterProDomain, ACMG_Predictors, sep = ":"),
                                        VSS2_classification = case_when(
                                          VSS2_score >= 10 ~ "Pathogenic",
                                          VSS2_score >= 5 & VSS2_score <= 9 ~ "Likely pathogenic",
                                          VSS2_score >= 0 & VSS2_score <= 4 ~ "VUS",
                                          VSS2_score < 0 ~ "Benign"),
                                        VSS2_version = "v2.3" #Change this to correpsonding version
                                      ))
  
  data_final_all <- future_map(data_final_all, ~ .x %>%
                                 mutate(
                                   ACMG = paste(ACMG_AF, ACMG_MutType, ACMG_Splicing, ACMG_InterProDomain, ACMG_Predictors, sep = ":"),
                                   VSS2_classification = case_when(
                                     VSS2_score >= 10 ~ "Pathogenic",
                                     VSS2_score >= 5 & VSS2_score <= 9 ~ "Likely pathogenic",
                                     VSS2_score >= 0 & VSS2_score <= 4 ~ "VUS",
                                     VSS2_score < 0 ~ "Benign"),
                                   VSS2_version = "v2.3" #Change this to correpsonding version
                                 ))
  
  
  ## VSS2_filter: Move ip columns to the front (apply to list) 
  col_index <- c(268,23,25,20,27,26,289,290,291,281,293,292,287,288,296,295,269,271,270,272,274,273,30,33,220,224,243,250,252,254,264,266,294)
  
  data_final_filtered_reorder <- lapply(data_final_filtered, function(df) {
    raw_cols <- colnames(df)
    
    valid_index <- col_index[col_index <= length(raw_cols)]
    cols_to_move <- raw_cols[valid_index]
    cols_rest <- raw_cols[!raw_cols %in% cols_to_move]
    
    # 重新排列：前12個保持不動，從第13個開始插入 cols_to_move
    new_order <- c(cols_rest[1:12], cols_to_move, cols_rest[-(1:12)])
    
    df[, new_order, drop = FALSE]
  })
  
  ## Anno_all: Move ip columns to the front
  data_final_all_reorder <- lapply(data_final_all, function(df) {
    raw_cols <- colnames(df)
    
    valid_index <- col_index[col_index <= length(raw_cols)]
    cols_to_move <- raw_cols[valid_index]
    cols_rest <- raw_cols[!raw_cols %in% cols_to_move]
    
    # 重新排列：前12個保持不動，從第13個開始插入 cols_to_move
    new_order <- c(cols_rest[1:12], cols_to_move, cols_rest[-(1:12)])
    
    df[, new_order, drop = FALSE]
  })
  
  
  ##Table 1: variant numbers by steps (1 sample 1 output file)
  data_final_filtered_rows <- as.data.frame(data_final_filtered_rows)
  
  for (i in 1:nrow(data_final_filtered_rows)) {
    file_name <- row.names(data_final_filtered_rows)[i]
    current_row <- data_final_filtered_rows[i, , drop = FALSE]
    save_name <- paste0("variant_numbers_by_steps_", file_name)
    
    write.table(current_row, file = file.path(output_path, save_name), row.names = F, col.names = T, quote = FALSE, sep = "\t")
  }
  
  ##Table 2: Anno_all (1 sample 1 output file)
  
  for (i in seq_along(data_final_all_reorder)) {
    
    df <- data_final_all_reorder[[i]]
    
    file_name <- names(data_final_all_reorder)[i]
    
    df_clean <- df %>%
      mutate(across(where(is.character), ~ gsub("[\t\r\n]", " ", .))) %>%
      mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse=","))))
    
    write.table(
      df_clean,
      file = file.path(output_path, paste0("Anno_all_", file_name)),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
    )
  }
  
  
  ##Table 3: VSS2_filter (1 sample 1 output file)
  for (i in seq_along(data_final_filtered_reorder)) {
    
    df <- data_final_filtered_reorder[[i]]
    
    file_name <- names(data_final_filtered_reorder)[i]
    
    df_clean <- df %>%
      mutate(across(where(is.character), ~ gsub("[\t\r\n]", " ", .))) %>%
      mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse=","))))
    
    write.table(
      df_clean,
      file = file.path(output_path, paste0("VSS2_filter_", file_name)),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
    )
  }
  
  ##### End of VSS2.3 main script #####
  
  # Save the entire workspace to an .RData file
  save_workspace()
  
}, error = function(e) {
  # Error handling block
  
  cat("Error occurred during script execution:\n")
  print(e)
  
  # Save the current state of R session
  save.image(file = file.path(output_path, "workspace_error_stage.RData"))
  cat("Current stage of workspace saved to", file.path(output_path, "workspace_error_stage.RData"), "\n")
})


