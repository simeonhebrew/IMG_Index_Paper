## IMG_Index_Paper

---
This repository contains scripts and source data useful to reproduce results obtained for the paper: "Immune surveillance as a governing principle of human gut microbiome structure"

---
### Description:

in `scripts`:

a) `IgA_binding.R` - an R script that contains code use to calculate the Kau index, probability score and probability ratio binding values for our cohorts as well as other public datasets that were integrated in the analysis.

b) `IM_score_calc.R` - an R script that contains code use to calculate the IMscore and IMG index for our cohorts which can be extended to other datasets.

c) `adm_method.R` - R script that contains the code that was used to determine the weights using the automatic democratic method involving inferring an optimal score for each species and performing non negative least square regression.

d) `bin_mags.sh`- a bash script that contains code for metagenomic assembly of sequencing reads using megahit followed by coverage estimation and binning using metabat2.

e) `metaphlan.sh` - a bash script that contains code for preprocessing of sequencing reads followed by taxonomic profiling using Metaphlan4.


in `source_data`:

Figure_2_SF, Figure_3_SF, Figure_4_SF, Figure_5_SF and Figure_6_SF contain the source data useful for reproducing the figure 2, 3, 4, 5 and 6 respectively.


### Software Installation:

a) To calculate the IgA binding, we advise installation of the [IgAscores](https://github.com/microbialman/IgAScores) R package authored by Jackson _et al_. 2021

b) To calculate the IM score and the IMG index, you can access our [IMGI](https://github.com/simeonhebrew/IMGMI) R package for all installation and usage details.



For any inquiry, feel free to open an issue or reach out to us on email at simeon.nthuku@inserm.fr or lejla.imamovic@inserm.fr
