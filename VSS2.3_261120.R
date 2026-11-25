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
###cmd: Rscript VSS2.2.R --data /path/to/your/vcf --ref /path/to/your/ref --output /path/to/results --BA1_restore --QUAL_filter

# Define the list of command-line options
option_list <- list(
  make_option(c("--data"), type = "character", help = "File path to the annotation files"),
  make_option(c("--ref"), type = "character", help = "File path to the reference tables"),
  make_option(c("--output"), type = "character", help = "File path to the output results"),
  make_option(c("--DP"), type = "double", default = 10, help = "Depth cutoff [default %default]", metavar = "Integer"),
  make_option(c("--AF"), type = "double", default = 0.2, help = "Allele frequency cutoff [default %default]", metavar = "Numerical"),
  make_option(c("--POP_AF"), type = "double", default = 0.05, help = "Rounded population allele frequency cutoff [default %default]", metavar = "Numerical"),
  make_option(c("--VSS2"), type = "double", default = 5, help = "VSS2 cutoff [default %default]", metavar = "Integer -2~13"),
  make_option(c("--QUAL_filter"), action = "store_true", default = FALSE, help = "Turn on QUAL filters (QD,SB,ReadPos) [default %default]"),
  make_option(c("--BA1_restore"), action = "store_true", default = FALSE, help = "Retain BA1 variants [default %default]"),
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
POP_AF_cutoff <- args$POP_AF
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
  cat("AF cutoff:", 0.2, "\n")
} else {
  cat("AF cutoff:", AF_cutoff, "\n")
}

if (is.null(POP_AF_cutoff)) {
  cat("Rounded POP_AF cutoff:", 0.05, "\n")
} else {
  cat("Rounded POP_AF cutoff:", POP_AF_cutoff, "\n")
}

if (is.null(VSS2_cutoff)) {
  cat("VSS2 cutoff:", 5, "\n")
} else {
  cat("VSS2 cutoff:", VSS2_cutoff, "\n")
}

if (!is.null(phenotype)) {
  cat("Phenotype filter:", paste(phenotype, collapse = ", "), "\n")
} else {
  cat("Phenotype filter: NULL\n")
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
  full.names = TRUE,
  ignore.case = TRUE
)

predictor_cal <- NULL
genes <- NULL
common <- NULL

for (f in csv_files) {
  
  if (grepl("Predictor_calibration", basename(f),
            ignore.case = TRUE)) {
    
    predictor_cal <- read.csv(f, check.names = FALSE)
    cat("predictor_calibration:", basename(f), "\n")
    
  } else if (grepl("Gene_panel_coordinates", basename(f),
                   ignore.case = TRUE)) {
    
    genes <- read.csv(f, check.names = FALSE)
    cat("gene_panels:", basename(f), "\n")
    
  } else if (grepl("Common_variants", basename(f),
                   ignore.case = TRUE)) {
    
    common <- read.csv(f, check.names = FALSE)
    cat("common_variants:", basename(f), "\n")
  }
}


### Check required and optional reference files


# Required
if (is.null(predictor_cal)) {
  stop("ERROR: Required file 'Predictor_calibration' not found.")
}

# Optional
if (is.null(genes)) {
  message("Optional file 'Gene_panel_coordinates' not found. Skipping.")
}

if (is.null(common)) {
  message("Optional file 'Common_variants' not found. Skipping.")
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
  
  ### ---------- 1. Load reference files----------
  
  csv_files <- list.files(
    ref_path,
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  # 統一讀取函數，避免重複寫 if/else
  load_reference <- function(prefix, required = FALSE) {
    
    files <- csv_files[
      grepl(paste0("^", prefix), basename(csv_files),
            ignore.case = TRUE)
    ]
    
    if (length(files) == 0) {
      if (required) stop("Required reference missing: ", prefix)
      message("Optional reference skipped: ", prefix)
      return(NULL)
    }
    
    # 避免同類型多個 CSV 被不小心覆蓋
    if (length(files) > 1) {
      stop("Multiple reference files found: ", prefix)
    }
    
    message("Loaded: ", basename(files))
    read.csv(files[1], check.names = FALSE)
  }
  
  predictor_cal <- load_reference("Predictor_calibration", TRUE)
  genes <- load_reference("Gene_panel_coordinates")
  common <- load_reference("Common_variants")
  
  
  ### ---------- 2. Clean raw data ----------

  # 整合 clean_raw_data + dedup_semicolon_cells
  clean_raw_data <- function(data) {
    
    data <- data %>%
      mutate(
        across(where(is.factor), as.character),
        across(where(is.character), ~ {
          
          s <- stri_replace_all_fixed(., "\\x3b", ";")
          s <- stri_replace_all_fixed(s, "\\x3d", "=")
          s <- stri_replace_all_fixed(s, ".,", "")
          
          parts <- stri_split_fixed(s, ";")
          
          vapply(parts, function(x) {
            
            if (all(is.na(x))) return(NA_character_)
            
            x <- stri_trim_both(x)
            
            # 同時移除缺失值及重複值
            x <- x[
              !is.na(x) &
                !(toupper(x) %in% c("", ".", "NA", "N/A", "NULL"))
            ]
            
            x <- unique(x)
            
            if (length(x) == 0) return(NA_character_)
            
            stri_join(x, collapse = ";")
            
          }, character(1))
        })
      )
    
    # OMIM 欄位不存在時自動略過
    if ("dbNSFP_OMIM_id" %in% names(data)) {
      
      data$dbNSFP_OMIM_id <- vapply(
        data$dbNSFP_OMIM_id,
        function(s) {
          
          if (is.na(s)) return(NA_character_)
          
          ids <- stri_trim_both(
            unlist(strsplit(s, "[,;]"))
          )
          
          # 僅格式化有效的 1–6 位數字
          if (!length(ids) ||
              !all(grepl("^[0-9]{1,6}$", ids))) {
            return(s)
          }
          
          paste(
            sprintf("%06d", as.integer(ids)),
            collapse = ";"
          )
        },
        character(1)
      )
    }
    
    return(data)
  }
  
  
  ###---------- 3. Merge equivalent annotation columns----------
  
  # 欄位不存在時不會報錯
  merge_columns <- function(data) {
    
    column_pairs <- c(
      verPhCons = "phastCons100way_vertebrate",
      verPhyloP = "phyloP100way_vertebrate",
      PHRED = "CADD_phred"
    )
    
    for (source in names(column_pairs)) {
      
      target <- unname(column_pairs[source])
      
      if (!source %in% names(data)) next
      
      if (!target %in% names(data)) {
        data[[target]] <- data[[source]]
      } else {
        idx <- is.na(data[[target]])
        data[[target]][idx] <- data[[source]][idx]
      }
      
      data[[source]] <- NULL
    }
    
    return(data)
  }
  
  
  ###---------- 4. Calculate genotype allele fraction----------

  add_gt_AF <- function(data, ad_col = "gt_AD") {
    
    # 檢查 AD 欄位
    if (!ad_col %in% names(data)) {
      stop("Missing AD column: ", ad_col)
    }
    
    data$gt_AF <- vapply(
      as.character(data[[ad_col]]),
      function(ad) {
        
        if (is.na(ad)) return(NA_real_)
        
        x <- suppressWarnings(
          as.numeric(strsplit(ad, ",", fixed = TRUE)[[1]])
        )
        
        # 僅計算 biallelic AD
        if (length(x) != 2 ||
            anyNA(x) ||
            any(x < 0) ||
            sum(x) == 0) {
          return(NA_real_)
        }
        
        x[2] / sum(x)
      },
      numeric(1)
    )
    
    return(data)
  }
  
  
  ####----------5. Assign gene panels (Optional) ----------

  check_panel_hits <- function(data, genes) {
    
    if (is.null(genes)) return(data)
    
    # 統一必要座標欄位名稱
    names(genes) <- sub(
      "^chrom$", "CHROM", names(genes), ignore.case = TRUE
    )
    names(genes) <- sub(
      "^start$", "Start", names(genes), ignore.case = TRUE
    )
    names(genes) <- sub(
      "^end$", "End", names(genes), ignore.case = TRUE
    )
    
    required <- c("CHROM", "Start", "End")
    
    if (!all(required %in% names(genes))) {
      stop("Gene panel missing CHROM, Start, or End.")
    }
    
    if (!all(c("CHROM", "POS") %in% names(data))) {
      stop("Variant data missing CHROM or POS.")
    }
    
    # 動態偵測所有非座標/Gene欄位
    exclude_cols <- c("CHROM", "START", "END", "GENE")
    
    panel_cols <- names(genes)[
      !(toupper(names(genes)) %in% exclude_cols)
    ]
    
    if (length(panel_cols) == 0) return(data)
    
    message("Detected panels: ",
            paste(panel_cols, collapse = ", "))
    
    # 統一染色體格式
    genes$CHROM <- sub(
      "^chr", "", as.character(genes$CHROM),
      ignore.case = TRUE
    )
    
    data_chr <- sub(
      "^chr", "", as.character(data$CHROM),
      ignore.case = TRUE
    )
    
    genes$Start <- as.numeric(as.character(genes$Start))
    genes$End <- as.numeric(as.character(genes$End))
    data_pos <- as.numeric(as.character(data$POS))
    
    genes[panel_cols] <- lapply(
      genes[panel_cols],
      function(x) suppressWarnings(
        as.numeric(as.character(x))
      )
    )
    
    # 檢查 Panel 值為 0/1/NA
    invalid <- vapply(
      genes[panel_cols],
      function(x) any(!is.na(x) & !(x %in% c(0, 1))),
      logical(1)
    )
    
    if (any(invalid)) {
      stop("Invalid panel values: ",
           paste(panel_cols[invalid], collapse = ", "))
    }
    
    # 所有 Panel 欄位預設為 0
    for (col in panel_cols) {
      data[[col]] <- 0L
    }
    
    # 先依染色體建立索引，減少重複搜尋
    chr_index <- split(
      seq_len(nrow(genes)),
      genes$CHROM
    )
    
    for (i in seq_len(nrow(data))) {
      
      if (is.na(data_chr[i]) || is.na(data_pos[i])) next
      
      idx <- chr_index[[data_chr[i]]]
      
      if (is.null(idx)) next
      
      hit_rows <- idx[
        which(
          genes$Start[idx] <= data_pos[i] &
            genes$End[idx] >= data_pos[i]
        )
      ]
      
      if (length(hit_rows) == 0) next
      
      for (col in panel_cols) {
        data[[col]][i] <- as.integer(
          any(genes[[col]][hit_rows] == 1, na.rm = TRUE)
        )
      }
    }
    
    return(data)
  }
  
  
  ###---------- 6. Assign common variants (Optional) ----------
  
  # 將 Common 獨立成 Function
  check_common_variants <- function(data, common) {
    
    if (is.null(common)) return(data)
    
    required <- c("CHROM", "POS", "REF", "ALT")
    
    if (!all(required %in% names(data)) ||
        !all(required %in% names(common))) {
      stop("Common variant annotation missing required columns.")
    }
    
    make_key <- function(df) {
      
      chr <- sub(
        "^chr", "", as.character(df$CHROM),
        ignore.case = TRUE
      )
      
      paste(chr, df$POS, df$REF, df$ALT, sep = "_")
    }
    
    # 排除不完整的 variant key
    valid_data <- complete.cases(data[, required])
    valid_common <- complete.cases(common[, required])
    
    data$Common <- as.integer(
      valid_data &
        make_key(data) %in% make_key(common[valid_common, ])
    )
    
    return(data)
  }
  
  
  ###---------- 7. Adding QD_filter, SB_filter, ReadPos_filter
  add_qual_filters <- function(df) {
    
    df %>%
      mutate(
        QD_filter = if_else(QD < 2, "Fail", "Pass"),
        SB_filter = if_else(FS > 60 | SOR > 3, "Fail", "Pass"),
        ReadPos_filter = if_else(
          is.na(ReadPosRankSum) | ReadPosRankSum >= -8,
          "Pass",
          "Fail"
        )
      )
  }
  
  
  ###---------- 8. Apply functions to all dataframes----------
  
  plan(multisession)
  
  table_list <- future_map(
    table_list,
    function(df) {
      
      # Data cleaning
      df <- clean_raw_data(df)
      df <- merge_columns(df)
      df <- add_gt_AF(df)
      
      # Optional annotations
      if (!is.null(genes)) {
        df <- check_panel_hits(df, genes)
      }
      
      if (!is.null(common)) {
        df <- check_common_variants(df, common)
      }
      
      # [新增] Quality filter annotation
      df <- add_qual_filters(df)
      
      return(df)
    },
    .options = furrr_options(
      packages = c("dplyr", "stringi")
    )
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
  ACMG_Splicing <- function(data, predictor_cal, sp_coltype = "Splicing_predictor") {
    
    # 取出 splicing predictor 欄位
    sp_cols <- predictor_cal %>% 
      filter(Coltype == sp_coltype) %>% 
      pull(Colname)
    
    # 定義單一值對照 cutoff 的 helper
    assign_splice_level <- function(value, colname) {
      if (is.na(value)) return(NA_character_)
      
      cutoffs <- predictor_cal %>%
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
  data3 <- future_map(data2, ~ACMG_Splicing(.x, predictor_cal, sp_coltype = "Splicing_predictor"))
  
  
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
  
  cols <- predictor_cal %>% filter(Coltype %in% c("Insilico_predictor", "Conservation_predictor")) %>% pull(Colname)
  
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
  
  score <- predictor_cal %>% filter(Coltype %in% c("Insilico_predictor", "Conservation_predictor"))
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

  
  ## Customized filter
  customized_filter <- function(
    df,
    DP_cutoff = 10,
    AF_cutoff = 0.2,
    POP_AF_cutoff = 0.05,
    VSS2_cutoff = 5,
    BA1 = TRUE,
    phenotype = NULL,
    QUAL_filter = FALSE
  ) {
    
    # 1. DP and AF filtering
    df <- df %>%
      filter(
        gt_DP >= .env$DP_cutoff,
        gt_AF >= .env$AF_cutoff
      )
    
    # 2. Population allele frequency filtering
    df <- df %>%
      filter(
        CLNID %in% .env$exclude_CLNID |
          if_all(all_of(af_cols), is.na) |
          if_any(
            all_of(af_cols),
            ~ round(.x, 2) <= .env$POP_AF_cutoff
          )
      )
    
    # 3. VSS2 score filtering
    df <- df %>%
      filter(VSS2_score >= .env$VSS2_cutoff)
    
    # 4. BA1 filtering
    if (!BA1) {
      df <- df %>%
        filter(is.na(ACMG_AF) | ACMG_AF != "BA1")
    }
    
    # 5. Quality filtering
    if (QUAL_filter) {
      df <- df %>%
        filter(
          QD_filter == "Pass",
          SB_filter == "Pass",
          ReadPos_filter == "Pass"
        )
    }
    
    # 6. Phenotype filtering
    if (!is.null(phenotype)) {
      df <- map_phenotype(df, phenotype)
    }
    
    return(df)
  }
  
  
  ## Apply customized filtering to all samples
  
  plan(multisession)
  
  data_final_filtered <- future_map(
    data_final_all,
    function(df) {
      
      customized_filter(
        df,
        DP_cutoff = DP_cutoff,
        AF_cutoff = AF_cutoff,
        POP_AF_cutoff = POP_AF_cutoff,
        VSS2_cutoff = VSS2_cutoff,
        BA1 = BA1_restore,
        QUAL_filter = QUAL_filter,
        phenotype = phenotype
      )
      
    },
    .options = furrr_options(
      packages = "dplyr"
    )
  )
  
  ### Count variants passing DP/AF and QUAL
  
  count_filter_steps <- function(df,
                                 DP_cutoff = 10,
                                 AF_cutoff = 0.2) {
    
    required_cols <- c(
      "gt_DP", "gt_AF",
      "QD_filter", "SB_filter", "ReadPos_filter"
    )
    
    if (!all(required_cols %in% names(df))) {
      stop("Missing columns: ",
           paste(setdiff(required_cols, names(df)),
                 collapse = ", "))
    }
    
    data.frame(
      DP10_AF0.2 = sum(
        df$gt_DP >= DP_cutoff &
          df$gt_AF >= AF_cutoff,
        na.rm = TRUE
      ),
      
      QUAL_pass = sum(
        df$QD_filter == "Pass" &
          df$SB_filter == "Pass" &
           df$ReadPos_filter == "Pass",
        na.rm = TRUE
      ))
  }
  
  
  ## Generate a report for the filtered variants
  
  # 從未經 customized filtering 的資料統計
  filter_counts <- lapply(
    data_final_all,
    count_filter_steps,
    DP_cutoff = 10,
    AF_cutoff = 0.2
  )
  
  data_final_filtered_rows <- data.frame(
    Original = sapply(table_list, nrow),
    DP10_AF0.2 = sapply(filter_counts, function(x) x$DP10_AF0.2),
    QUAL_pass = sapply(filter_counts, function(x) x$QUAL_pass),
    ACMG_AF = sapply(data1, function(df) sum(!is.na(df$ACMG_AF))),
    ACMG_MutType = sapply(data2, function(df) sum(!is.na(df$ACMG_MutType))),
    ACMG_Splicing = sapply(data3, function(df) sum(!is.na(df$ACMG_Splicing))),
    ACMG_InterProDomain = sapply(data4, function(df) sum(!is.na(df$ACMG_InterProDomain))),
    ACMG_Predictors = sapply(data5, function(df) sum(!is.na(df$ACMG_Predictors))),
    Variants_Filtered = sapply(data_final_filtered, nrow),
    check.names = FALSE
  )
  
  print(data_final_filtered_rows)
  
  #####----------Step 8. Output result files----------
  
  # Add ACMG combination & VSS2_classification/version
  data_final_filtered <- future_map(data_final_filtered, ~ .x %>%
                                      mutate(
                                        ACMG = paste(ACMG_AF, ACMG_MutType, ACMG_Splicing, ACMG_InterProDomain, ACMG_Predictors, sep = ":"),
                                        VSS2_classification = case_when(
                                          VSS2_score >= 10 ~ "Pathogenic",
                                          VSS2_score >= 6 & VSS2_score <= 9 ~ "Likely pathogenic",
                                          VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
                                          VSS2_score <= 1 ~ "Benign"),
                                        VSS2_version = "v2.3" #Change this to correpsonding version
                                      ))
  
  data_final_all <- future_map(data_final_all, ~ .x %>%
                                 mutate(
                                   ACMG = paste(ACMG_AF, ACMG_MutType, ACMG_Splicing, ACMG_InterProDomain, ACMG_Predictors, sep = ":"),
                                   VSS2_classification = case_when(
                                     VSS2_score >= 10 ~ "Pathogenic",
                                     VSS2_score >= 6 & VSS2_score <= 9 ~ "Likely pathogenic",
                                     VSS2_score >= 2 & VSS2_score <= 5 ~ "VUS",
                                     VSS2_score <= 1 ~ "Benign"),
                                   VSS2_version = "v2.3" #Change this to correpsonding version
                                 ))
  
  
  ## VSS2_filter: Move ip columns to the front (apply to list) 
  col_index <- c(268,13,23,25,20,27,26,274,275,276,269,270,271,273,272,289,290,291,283,292,294,295,30,33,220,224,243,250,252,254,264,266,293)
  
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
  
  
  ## Output 1: Variant numbers by steps (1 sample 1 output file)
  data_final_filtered_rows <- as.data.frame(data_final_filtered_rows)
  
  for (i in 1:nrow(data_final_filtered_rows)) {
    file_name <- row.names(data_final_filtered_rows)[i]
    current_row <- data_final_filtered_rows[i, , drop = FALSE]
    save_name <- paste0("variant_numbers_by_steps_", file_name, ".table")
    
    write.table(current_row, file = file.path(output_path, save_name), row.names = F, col.names = T, quote = FALSE, sep = "\t")
  }
  
  ## Output 2: Anno_all (1 sample 1 output file)
  
  for (i in seq_along(data_final_all_reorder)) {
    
    df <- data_final_all_reorder[[i]]
    
    file_name <- names(data_final_all_reorder)[i]
    
    df_clean <- df %>%
      mutate(across(where(is.character), ~ gsub("[\t\r\n]", " ", .))) %>%
      mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse=","))))
    
    write.table(
      df_clean,
      file = file.path(output_path, paste0("Anno_all_", file_name, ".table")),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
    )
  }
  
  
  ## Output 3: VSS2_filter (1 sample 1 output file)
  for (i in seq_along(data_final_filtered_reorder)) {
    
    df <- data_final_filtered_reorder[[i]]
    
    file_name <- names(data_final_filtered_reorder)[i]
    
    df_clean <- df %>%
      mutate(across(where(is.character), ~ gsub("[\t\r\n]", " ", .))) %>%
      mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse=","))))
    
    write.table(
      df_clean,
      file = file.path(output_path, paste0("VSS2_filter_", file_name, ".table")),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
    )
  }
  
  ##### End of VSS2.2 main script #####
  
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


