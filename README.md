# Gene Set Enrichment Analysis (GSEA) Script

**Author:** Jonas D. Elsborg  
**Email:** [jonas.elsborg@path.ox.ac.uk](mailto:jonas.elsborg@path.ox.ac.uk)  

---

## Overview

This R script performs **gene set enrichment analysis (GSEA)** using Fisher's exact test on categorical data such as Gene Ontology (GO) annotations.  
It is designed to work with proteomics or other datasets containing annotated protein groups and a logical column defining a *foreground* (e.g., significantly regulated proteins).

The script:
- Parses GO annotation columns containing multiple terms per entry (semicolon-separated).
- Computes enrichment statistics (odds ratio, enrichment factor).
- Optionally applies the **modified Haldane-Anscombe correction** to handle zero counts and infinite odds ratios.
- Adjusts p-values using the **Benjamini-Hochberg (FDR)** correction.
- Returns a tidy dataframe summarizing enrichment per term and GO category (BP, MF, CC).

---

## Dependencies

The following R packages are required:

```r
tidyverse
```

## Input Data

The script expects a tab-separated file:

protein_groups_annotated.txt

This file should contain:

One or more columns with GO terms (e.g., Gene Ontology (biological process), Gene Ontology (molecular function), etc.).

A column with UniProt IDs (semicolon-separated if multiple).

A logical column (e.g. foreground) marking foreground vs background proteins.

Example snippet:

|  Protein IDs  | Gene Ontology (biological process) | Gene Ontology (molecular function) | foreground |
|:-------------:|:----------------------------------:|------------------------------------|------------|
| P12345;Q67890 | DNA repair; Chromatin organization | ATP binding                        | TRUE       |
| Q99999        | RNA processing                     | Nucleotide binding                 | FALSE      |

## Usage

Place the script and protein_groups_annotated.txt in the same folder.

Open the script in RStudio and run it - it automatically sets the working directory to the script's location.

The output file will be written to the same directory.

## Function: go_enrichment()
```
go_enrichment(
  df,               # Dataframe with GO terms and a foreground column
  category_cols,    # Columns containing GO annotations (e.g., starts_with("Gene Ontology"))
  foreground_col,   # Boolean column marking foreground entries
  sep = ";\\s*",    # Separator for GO terms, default is ;
  mHA_corr = FALSE  # Apply modified Haldane-Anscombe correction.
)
```
### Output columns

| Column | Description |
|---------|--------------|
| `GO_type` | GO category (BP, MF, CC) |
| `GO_term` | The individual GO term |
| `a`, `b`, `c`, `d` | Counts from the contingency table |
| `odds_ratio` | Fisher’s exact test odds ratio |
| `enrichment_factor` | Ratio of observed vs expected probability |
| `direction` | “enriched” or “depleted” |
| `p_value` | Raw Fisher test p-value |
| `p_adj` | FDR-adjusted p-value |

## Output

enrichment_fisher_test.txt

Tab-separated results of the enrichment analysis.

## Reference

Modified Haldane-Anscombe correction described in:

F. Weber, G. Knapp, K. Ickstadt, G. Kundt, Ä. Glass, *Zero-cell corrections in random-effects meta-analyses.* Res. Synth. Methods 11, 913-919 (2020).
  

## Example Call
```
enrichment_result <- go_enrichment(
  go_annotated,
  category_cols = starts_with("Gene Ontology"),
  foreground_col = foreground,
  mHA_corr = TRUE
)

enrichment_result %>% 
  write_tsv("enrichment_fisher_test.txt")
```
## Notes

It is recommended to enable the mHA_corr = TRUE option to avoid infinite odds ratios.

The script abbreviates long GO column names automatically (BP, MF, CC).

Works with any categorical annotation columns, not limited to GO terms.

### License

This script is released under the MIT License.
You are free to use and modify it with proper attribution.
