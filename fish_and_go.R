# Script perform geneset enrichment analysis
#
# By Jonas D. Elsborg
# jonas.elsborg@path.ox.ac.uk
# Oct25

# Description:


# clear env
rm(list = ls())

# set wd
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)

# import dependencies
library(tidyverse)
#library(hablar)

#### read data
# read dataframe, protein groups with uniprot ID's in a column. 
# The function will take the first of several uniprot ID's per row, separated by semi-colon. 
# If you prefer, you can choose the uniprot ID that you think best matches your protein of interest.

go_annotated <- read_tsv(file = "protein_groups_annotated.txt")

### GSEA function
# It handles infinite odds ratios using the modified Haldane-Anscombe correction if this is desired
# https://doi.org/10.1002/jrsm.1460
# 
# It will perform GSEA on any catagorical columns (they should be supplied as factors). 
# If they are GO-terms, the names will be abbreviated.
# P-values are adjusted using the Benjamini-Hochberg correction
# Besides odds ratio, an enrichment factor is calculated which reports ratio for the expected probability under random sampling from the full dataset.
# The function also needs a true/false column to indicate the foreground from the background.

go_enrichment <- function(df, category_cols, foreground_col, sep = ";\\s*", mHA_corr = FALSE) {
  foreground_col <- rlang::ensym(foreground_col) # capture column name (tidy eval)
  ha <- 0.5 # Haldane-Anscombe correction factor
  
  # 1) Pivot GO columns longer and split semicolon-separated terms
  df_long <- df %>%
    pivot_longer(
      cols = {{ category_cols }},  # accepts tidyselect (e.g. starts_with("Gene Ontology"))
      names_to = "GO_type", # names (BF/CC/MF or other)
      values_to = "GO_term" # values (the GO-terms)
    ) %>%
    filter(!is.na(GO_term) & GO_term != "") %>% # purge empty rows
    separate_rows(GO_term, sep = sep) %>% # split GO-terms into multiple rows
    mutate(GO_type = case_when( # renames long GO names
      str_detect(GO_type, "biological") ~ "BP",
      str_detect(GO_type, "molecular") ~ "MF",
      str_detect(GO_type, "cellular") ~ "CC",
      TRUE ~ GO_type
    ))
  
  # 2) Compute counts and run Fisher's exact test per GO term
  enrichment <- df_long %>%
    group_by(GO_type, GO_term) %>% # each term
    summarise(
      a = sum(!!foreground_col),              # foreground & in term
      c = sum(! (!!foreground_col)),          # background & in term
      .groups = "drop"
    ) %>%
    # b and d are calculated from total foreground/background sizes (taken from the original df)
    mutate(
      b = sum(df[[rlang::as_string(foreground_col)]]) - a,  # foreground & not in term
      d = sum(!df[[rlang::as_string(foreground_col)]]) - c  # background & not in term
    ) %>%
    # Rowwise so we can run fisher.test one GO term at a time and store the HT test object
    rowwise() %>%
    mutate( 
      # Run Fisher's exact test for the 2x2 table: [ [a b], [c d] ]
      test = list(fisher.test(matrix(c(a, b, c, d), nrow = 2))),
      p_value = test$p.value,
      odds_ratio =  test$estimate, # ((a * d) / (b * c)), # with numbers it is the true odds ratio. with test it is Conditional MLE (hypergeometric-based)
      enrichment_factor = (a / (a + b)) / ((a + c) / (a + b + c + d)), # enrichment factor
      # labels for plotting
      direction = case_when(
        is.na(odds_ratio) ~ NA_character_,
        odds_ratio > 1 ~ "enriched",
        odds_ratio < 1 ~ "depleted",
        TRUE ~ "neutral"
      )
    ) %>%
    ungroup() %>%
    # 3) Correct p-values per GO_type (BH correction within BP, MF, CC separately)
    group_by(GO_type) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(GO_type, p_adj) %>%
    # 4) clean up ordering and select the useful columns
    select(GO_type, GO_term, a, b, c, d, odds_ratio, enrichment_factor,
           direction, p_value, p_adj)
  
  # 5) apply “Haldane-Anscombe” (HA) correction at rows with any 0
  if(mHA_corr){
    enrichment <- enrichment %>%
      rowwise() %>%
      mutate(
        # Check if any count is zero
        zero_flag = any(c(a, b, c, d) == 0),
        
        # Apply HA correction if needed
        a_ha = ifelse(zero_flag, a + ha, a),
        b_ha = ifelse(zero_flag, b + ha, b),
        c_ha = ifelse(zero_flag, c + ha, c),
        d_ha = ifelse(zero_flag, d + ha, d),
        
        # Recalculate odds ratio and enrichment factor
        odds_ratio = (a_ha * d_ha) / (b_ha * c_ha),
        enrichment_factor = (a_ha / (a_ha + b_ha)) / ((a_ha + c_ha) / (a_ha + b_ha + c_ha + d_ha))
      ) %>%
      ungroup()
  }
  
  
  enrichment # the df object within the function 
}


enrichment_result <- go_enrichment(
  go_annotated, # The dataframe, with GO terms
  category_cols = starts_with("Gene Ontology"), # This takes several columns, or you just reference a column directly 
  foreground_col = foreground, # Boolean column, TRUE denotes foreground
  mHA_corr = TRUE # It is recommended to enable this function to apply the modified 
)

enrichment_result %>% write_tsv(file = "enrichment_fisher_test.txt")