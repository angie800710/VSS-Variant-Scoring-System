# VSS: Variant Scoring System
VSS is an integrated Variant Scoring System for automated pathogenic variant prioritization. VSS scores ranging from −2 to 13, with higher scores indicating greater pathogenic potential. Score cutoffs of 3–6 represented the most favorable operating range for prioritizing variants with pathogenic potential. For categorical interpretation, scores ≥6 are recommended to correspond to Pathogenic/Likely Pathogenic (P/LP), scores of 2–5 to Variants of Uncertain Significance (VUS), and scores ≤1 to Benign/Likely Benign (B/LB).

## Annotation resources

VSS integrates multiple annotation sources through ANNOVAR. Users should first download and install ANNOVAR from:
https://www.openbioinformatics.org/annovar/annovar_download_form.php

Because the complete set of annotation databases used by VSS can require up to approximately 1 TB of storage, these resources are not hosted directly in this GitHub repository. Users who wish to access the curated annotation resources and configurations used by VSS may contact em61330@ncku.edu.tw to obtain the download link.

We recommend using the curated VSS annotation resources before running the VSS pipeline, as the workflow expects specific annotation fields and column names. Differences in database versions, annotation configurations, or column naming may affect compatibility with the pipeline.

**VSS currently uses the following annotation databases and resources**
- **CADD** — Combined Annotation Dependent Depletion, v1.7 (February 2024)
- **dbscSNV** — Database of Splicing Consensus Single-Nucleotide Variants, v1.1 (December 2015)
- **gnomAD** — Genome Aggregation Database, v4.1 (March 2024)
- **RefSeq Gene (refGene)** — RefSeq gene annotation (August 2022)
- **ClinVar** — Database of clinical significance of variants (June 2025)
- **dbNSFP** — Functional prediction and annotation database for nonsynonymous single-nucleotide variants, v5.1a (March 2025)
- **Ensembl Gene (ensGene)** — Ensembl gene annotation based on GENCODE v46 (October 2024)
- **SpliceAI** — Deep learning–based splice variant prediction resource, v1.3 (October 2019)
- **Taiwanese MAF** — Minor allele frequency database derived from the Taiwan Biobank (July 2023)
- **dbSNP** — NCBI database of genetic variation, build 156 (July 2023)

## VSS implementation
VSS is implemented in R and executed through an automated command-line interface. The input is the ANNOVAR-annotated vcf, and the VSS workflow is executed using the following command: ```
Rscript VSS.R --data <annotated_vcf> --ref <reference_directory> --output <output_directory>
```. The reference directory is provided in this repository and contains the required ```Used_columns_scores.csv``` file for in-silico predictor standardization. Two additional reference files, ```Gene_panel_coordinates.csv``` and ```Common_variants_after_mask.csv```, are optional and may be used for user-defined gene panel filtering and variant masking, respectively.

## VSS parameters
The following parameters are available for user-defined filtering:

| Parameter | Description | Default |
|---|---|---|
| --DP | Minimum sequencing depth (DP) required for variant retention | 10 |
| --AF | Minor allele frequency (MAF) threshold | 0.05 |
| --VSS | VSS score threshold for variant prioritization | 6 |
| --QUAL_filter | Enable additional quality filters based on quality depth, strand bias, and read-position bias | Enabled |
| --BA1_restore | Restore variants classified as BA1 (stand-alone benign evidence) | Disabled |
| --phenotype | Phenotype keyword filtering using comma-separated terms.* | None |

\* VSS supports partial matching of phenotype descriptions specified using the `--phenotype` option. Multiple phenotype terms are separated by commas, with underscores used in place of spaces within each term. Matched HPO terms are reported in `Phenotype_match`, and the number of matched terms is recorded in `Phenotype_hit`.

## VSS outputs
The VSS pipeline generates three output tables for each analyzed sample, together with an R workspace file:

| Output file | Description |
|---|---|
| `Anno_all_sampleID` | Contains the complete set of annotated variants after integration of ANNOVAR annotations and VSS-specific evidence features. |
| `VSS_filter_sampleID` | Contains variants retained after applying VSS filtering according to the specified command-line parameters. |
| `variant_numbers_by_steps_sampleID` | Summarizes the number of variants processed at each VSS module. |
| `R workspace` | Preserves the complete analysis environment, including intermediate objects, parameter settings, execution records, and error messages. |


## Dependencies
VSS is implemented in **R** and requires the following R packages:
- `dplyr`
- `tidyr`
- `magrittr`
- `furrr`
- `purrr`
- `stringr`
- `reshape2`
- `optparse`
- `stringi`
- `data.table`
- `vcfR`

The required packages can be installed in R using:

```r
install.packages(c("dplyr","tidyr","magrittr","furrr","purrr","stringr","reshape2","optparse","stringi","data.table","vcfR"))
```

## Publication
For detailed information, please refer to the following paper.
VSS: An Integrated Variant Scoring System for Automated Pathogenic Variant Prioritization.
Ying-An Chen, Yen-An Tang, Chiao-May Chang, and H. Sunny Sun


