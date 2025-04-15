# Genome-Wide Association and Genomic Prediction in Yellow Bean Collection (YBC)

This repository contains code and data for the publication:
**"Genome-wide association and genomic prediction for Fe and Zn concentration and Fe bioavailability in a yellow bean collection of dry beans"**  
🔗 [Frontiers in Genetics, DOI: 10.3389/fgene.2024.1330361](https://doi.org/10.3389/fgene.2024.1330361)

---

## 🧬 Abstract

Dry bean (*Phaseolus vulgaris* L.) is a nutrient-dense crop central to global diets and the focus of biofortification efforts. This study explores the genetic basis and predictive modeling of iron (Fe), zinc (Zn), and Fe bioavailability (FeBio) in a 295-genotype **Yellow Bean Collection (YBC)** evaluated across two years and locations.

Key insights:
- GWAS revealed trait-specific quantitative trait nucleotides (QTNs).
- Genomic prediction (GP) accuracies ranged from 0.12 (Fe) to 0.72 (FeBio).
- Adding GWAS-derived QTNs improved FeBio prediction by 0.03.
- Fe concentration and FeBio were largely uncorrelated, supporting FeBio as an independent trait for breeding.

---

## 🗂 Repository Structure

### 📁 Scripts

| File | Description |
|------|-------------|
| `1.Phenotype_varianCompo_GWAS_popStruct_YBC.R` | Analyzes correlations, performs variance component analysis, PCA, kinship matrix calculation, GWAS, and visualizations for phenotypes and seed color |
| `2.RKHS_SSI_Yield_Fe_Zn_YBC.R` | Runs single-site inference (SSI) and RKHS models for Yield, Fe, and Zn across MI and NE (2018–2019) |
| `3.RKHS_SSI_FeBio_YBC.R` | Fits RKHS models for FeBio including GWAS-assisted fixed effects |
| `4.RKHS_SSI_YBC_predictionSet.R` | Predicts trait values in the full Yellow Bean Collection using RKHS/SSI models |
| `5.RKHS_SSI_Andean_predictionSet.R` | Applies models specifically to Andean accessions in the YBC prediction set |

---

### 📁 Data

- `YBC_GWASAssistedGP.RData`:  
  Contains phenotype, genotype, hapmap, and SNP position data for all analyses.

- `YBC_phenotype_Color_CT.csv`:  
  Contains agronomic, mineral, and seed quality traits, cooking time, and population structure annotations.  
  Metadata was curated from:
  - [Sadohara et al., 2021](https://doi.org/10.1002/tpg2.20173)
  - [Sadohara et al., 2022](https://doi.org/10.1007/s10722-021-01323-0)

- `GBS_barcode_plate_info_YBC.txt`:  
  Plate layout and metadata for genotyping-by-sequencing samples.

🔗 **Raw sequencing data**:  
Available through [NCBI BioProject PRJNA1061170](https://www.ncbi.nlm.nih.gov/bioproject/1061170)

---

## 📦 Software Requirements

- **R** (≥ 4.0.0)
- Key packages:
  - `tidyverse`
  - `data.table`
  - `rrBLUP`
  - `sommer`
  - `BGLR`
  - `FactoMineR`, `ggplot2` (for PCA, plots)

---

## ✅ Outputs

- Trait variance components
- GWAS summary tables and Manhattan plots
- Prediction accuracy metrics
- Marker effects and model coefficients

---

## 📣 Citation

If you use this repository, please cite:

> Izquierdo, P., et al. (2024).  
> *Genome-wide association and genomic prediction for Fe and Zn concentration and Fe bioavailability in a yellow bean collection of dry beans.*  
> Frontiers in Genetics. [DOI: 10.3389/fgene.2024.1330361](https://doi.org/10.3389/fgene.2024.1330361)

---
