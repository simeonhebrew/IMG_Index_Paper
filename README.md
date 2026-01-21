## IMG_Index_Paper

---
This repository contains scripts used to reproduce results obtained for the paper: "Immune surveillance as a governing principle of human gut microbiome structure"

---
Description:

in `scripts`:
a) `IgA_binding.R` - an R script that contains code use to calculate the Kau index, probability score and probability ratio binding values for our cohorts as well as other public datasets that were integrated in the analysis.

b) `adm_method.R` - R script that contains the code that was used to determine the weights using the automatic democratic method involving inferring an optimal score for each species and performing non negative least square regression.

c) `bin_mags.sh`- a SLURM script that contains code for metagenomic assembly of sequencing reads using megahit followed by coverage estimation and binning using metabat2.

d) `metaphlan.sh` - a SLURM script that contains code for preprocessing of sequencing reads followed by taxonomic profiling using Metaphlan4.


