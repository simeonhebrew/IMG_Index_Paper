#To calculate the IM score and IMG index, install the IMG index R package as detailed in the README file.
library(IMGI)
library(IgAScores)


#The user requires a table of taxa against samples that can be obtained from phyloseq or any other general 
microbiome processing

#Example for HC and CIS IMG index calculation


#Calculating IM score and global score in HC
otu_pos_hc_df <- as.data.frame(otu_table(phyloseq_hc_pos))

rownames(otu_pos_hc_df) <- probratio_scores_hc_pos_df$Genus

otu_pos_hc_df_rel <- relabund(otu_pos_hc_df)

im_score_hc_pos <- calculate_im_score(otu_pos_hc_df_rel, probratio_scores_hc_pos)

global_score_hc_pos <- calculate_im_global_score(im_score_hc_pos)

rownames(global_score_hc_pos) <- paste0("HC_", rownames(global_score_hc_pos))



#Calculating IM score and global score in CIS
otu_pos_cis_df <- as.data.frame(otu_table(phyloseq_cis_pos))

rownames(otu_pos_cis_df) <- probratio_scores_cis_pos_df$Species

otu_pos_cis_df_rel <- relabund(otu_pos_cis_df)

im_score_cis_pos <- calculate_im_score(otu_pos_cis_df_rel, probratio_scores_cis_pos)

global_score_cis_pos <- calculate_im_global_score(im_score_cis_pos)

rownames(global_score_cis_pos) <- paste0("CIS_", rownames(global_score_cis_pos))



#Calculating IMG index 

global_score_hc_cis_pos <- calculate_isgm_index(global_score_hc_pos, global_score_cis_pos)


