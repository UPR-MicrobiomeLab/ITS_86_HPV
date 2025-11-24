# README

This dataset includes:
- A curated **sample file** (metadata) 
- A processed **phyloseq object** derived from ITS data  
- R scripts to reproduce the full analysis pipeline  

---

## Phyloseq Object

The processed phyloseq object is included:

**File:** `PHYraw_SH_ITS_alpha.RData`  

---

## Scripts

We included the R scripts used to generate the **phyloseq objects** from raw QIIME files, as well as to perform the **alpha** and **beta diversity analyses**.

These scripts contain all commands and workflow used to generate:

- Alpha diversity metrics (Shannon, Simpson, Observed Richness)  
- Beta diversity dissimilarities (Bray–Curtis, Jaccard, Aitchison CLR)  
- Filtering steps  
- Metadata merging  
- Statistical tests (mixed models, PERMANOVA, etc.)  
- Final exported tables and figures  

**Available files:**
- `1_generate_phyloseq_obj_ITS.R`
- `alpha_diversity_script.R`
- `beta_diversity_script.R`

These scripts can be used as **templates** to reproduce the analyses.

---

