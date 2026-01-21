#Inferring IgA binding of species in all represented cohorts
library(IgAScores)

#HC


rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_hc_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_hc_neg))



rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Species, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Species, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- tax_pos_samples$Species
rownames(otu_neg_samples_df_inv) <- tax_neg_samples$Species



colnames(otu_pos_samples_df_inv) <- gsub("_IgApos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_IgAneg", "", colnames(otu_neg_samples_df_inv))



rel_otu_pos <- relabund(otu_pos_samples_df_inv)
rel_otu_neg <- relabund(otu_neg_samples_df_inv)




#Fraction sizes 
possize_control <- c(CER005 = 0.091,CER203=0.0319, CER207=0.0245, CER209=0.045, 
CER210=0.0462, CER211=0.1162, CER212=0.0123, CER215=0.028, CER216=0.037, CER218=0.059, 
CER220=0.0276, CER223=0.0327, CER224=0.0392, CER226=0.0444, CER231=0.0501,CER232=0.0414, 
CER233=0.053, CER235=0.0609)

negsize_control <- c(CER005 = 0.3746, CER203=0.3038, CER207=0.191, CER209=0.3528, 
CER210=0.3161, CER211=0.3829, CER212=0.4059, CER215 = 0.201, CER216=0.313, CER218=0.2832, 
CER220=0.4229, CER223=0.2659, CER224=0.3547, CER226=0.3806, CER231=0.4411, CER232=0.2932, 
CER233=0.4073, CER235=0.2539)



lowest <- min(min(rel_otu_pos[rel_otu_pos!=0]),min(rel_otu_neg[rel_otu_neg!=0]))
lowest

pseudo = 9.006291e-07





prob_ratio_pos_neg_control_filtered <- igascores(posabunds = rel_otu_pos, negabunds = 
rel_otu_neg,
                                                 possizes = possize_control, negsizes = 
negsize_control,
                                                 pseudo=pseudo)



#CIS

otu_pos_samples_df <- as.data.frame(t(otu_table(phyloseq_cis_pos)))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(t(otu_table(phyloseq_cis_neg)))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)




rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_cis_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_cis_neg))



rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Species, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Species, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- tax_pos_samples$Species
rownames(otu_neg_samples_df_inv) <- tax_neg_samples$Species



colnames(otu_pos_samples_df_inv) <- gsub("_IgApos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_IgAneg", "", colnames(otu_neg_samples_df_inv))



rel_otu_pos <- relabund(otu_pos_samples_df_inv)
rel_otu_neg <- relabund(otu_neg_samples_df_inv)



possize_cis <- c(ENT020=0.0098, ENT042=0.0402, ENT044=0.0136, ENT049=0.044, 
ENT050=0.0077,ENT073=0.0387,ENT079=0.0536,ENT101=0.0393,ENT103=0.0857, ENT112=0.0523, 
ENT113=0.041, ENT115=0.0728, ENT116 = 0.0156, ENT117=0.0634)

negsize_cis <- c(ENT020=0.4304, ENT042=0.323, ENT044=0.2605, ENT049=0.4101, 
ENT050=0.3041,ENT073=0.2986,ENT079=0.7219,ENT101=0.5114,ENT103=0.4002, ENT112=0.3034, 
ENT113=0.4297, ENT115=0.3551, ENT116=0.5825, ENT117=0.4097)




lowest <- min(min(rel_otu_pos[rel_otu_pos!=0]),min(rel_otu_neg[rel_otu_neg!=0]))
lowest

pseudo = 2.708015e-06





prob_ratio_pos_neg_cis_filtered <- igascores(posabunds = rel_otu_pos, negabunds = 
rel_otu_neg,
                                             possizes = possize_cis, negsizes = 
negsize_cis,
                                             pseudo=pseudo)




#Olm et al.


otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_olm_sp_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_olm_sp_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_olm_sp_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_olm_sp_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Species, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Species, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Family, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Family, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Phylum, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Phylum, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Class, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Class, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("_p", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_r", "", colnames(otu_neg_samples_df_inv))


otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]



otu_neg_samples_df_inv$"2033" <- NULL
otu_pos_samples_df_inv$"2033" <- NULL

rel_otu_pos <- relabund(otu_pos_samples_df_inv)
rel_otu_neg <- relabund(otu_neg_samples_df_inv)


colnames(rel_otu_pos) <- gsub("X", "", colnames(rel_otu_pos))
colnames(rel_otu_neg) <- gsub("X", "", colnames(rel_otu_neg))



pos_size =c("1022"=0.466, "1026"=0.323, "1036"=0.563, "1037"=0.545, "1038"=0.751, 
"1041"=0.684, "2023"=0.488, "2024"=0.566, "2025"=0.526, "2027"=0.703, "2028"=0.681,
            "2029"=0.512, "2030"=0.682, "2034"=0.526, "2039"=0.719, "2040"=0.471, 
"7022"=0.417, "7023"=0.414, "7024"=0.531, "7028"=0.807, "7029"=0.553,
            "7033"=0.601, "7035"=0.309, "7037"=0.790, "7039"=0.549, "7040"=0.740, 
"7041"=0.431, "8018"=0.739, "8026"=0.645, "8030"=0.580, "8032"=0.602, "8038"=0.676)


neg_size =c("1022"=0.1353, "1026"=0.0670, "1036"=0.0726, "1037"=0.1414, "1038"=0.0398, 
"1041"=0.0498, "2023"=0.0286, "2024"=0.1131, "2025"=0.0085, "2027"=0.0268, "2028"=0.0954,
            "2029"=0.0948, "2030"=0.0751, "2034"=0.0947, "2039"=0.0234, "2040"=0.0588, 
"7022"=0.1747, "7023"=0.0949, "7024"=0.1427, "7028"=0.0767, "7029"=0.1251,
            "7033"=0.0627, "7035"=0.0290, "7037"=0.0800, "7039"=0.0855, "7040"=0.1344, 
"7041"=0.0372, "8018"=0.0961, "8026"=0.1110, "8030"=0.0627, "8032"=0.1080, "8038"=0.0860)




lowest <- min(min(rel_otu_pos[rel_otu_pos!=0]),min(rel_otu_neg[rel_otu_neg!=0]))
lowest

pseudo = 1.000807e-07

#Iga probability ratio score 

prob_scores_pos_neg_olm <- igascores(posabunds = rel_otu_pos, negabunds = rel_otu_neg,
                                     possizes = pos_size, negsizes = neg_size,
                                     pseudo=pseudo)




#van Gogh et al.

otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_merel_sp_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_merel_sp_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_merel_sp_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_merel_sp_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Species, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Species, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Family, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Family, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Phylum, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Phylum, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Class, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Class, unique = TRUE)


otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]



colnames(otu_pos_samples_df_inv) <- gsub("_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_pre", "", colnames(otu_neg_samples_df_inv))



rel_otu_pos <- relabund(otu_pos_samples_df_inv)
rel_otu_neg <- relabund(otu_neg_samples_df_inv)


colnames(rel_otu_pos) <- gsub("X", "", colnames(rel_otu_pos))
colnames(rel_otu_neg) <- gsub("X", "", colnames(rel_otu_neg))





pos_size =c("001"=0.6665, "002"=0.29, "003"=0.3355, "004"=0.2145, "005"=0.805, "006"=0.648, 
"007"=0.5885,
            "008"=0.44, "009"=0.4015, "010"=0.481)




lowest <- min(min(rel_otu_pos[rel_otu_pos!=0]),min(rel_otu_neg[rel_otu_neg!=0]))
lowest

#pseudo = 2.001351e-07


#Iga probability score
prob_only_scores_pos_neg_merel <- igascores(posabunds = rel_otu_pos, possizes =  pos_size, 
presortabunds = rel_otu_neg, method = "prob")



######Probstel et al. Remitting-Relapsing Multiple sclerosis cohort############
#HC

otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_hc_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_hc_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_hc_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_hc_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("_IgA_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_IgA_neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_hc <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_hc <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_hc[rel_otu_pos_hc!=0]),min(rel_otu_neg_hc[rel_otu_neg_hc!=0]))
lowest


#pseudo =  4.193294e-06



kau_scores_probstel_hc<- igascores(posabunds = rel_otu_pos_hc, negabunds = rel_otu_neg_hc, 
pseudo =  4.193294e-06,  method = "kau")


#Remission
otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_rem_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_rem_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_rem_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_rem_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("_IgA_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_IgA_neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_rem <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_rem <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_rem[rel_otu_pos_rem!=0]),min(rel_otu_neg_rem[rel_otu_neg_rem!=0]))
lowest


#pseudo = 5.700052e-06



kau_scores_probstel_rem<- igascores(posabunds = rel_otu_pos_rem, negabunds = 
rel_otu_neg_rem, pseudo =  5.700052e-06,  method = "kau")



#Relapse
otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_rel_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_probstel_ms_rel_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_rel_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_probstel_ms_rel_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("_IgA_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_IgA_neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_rel <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_rel <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_rel[rel_otu_pos_rel!=0]),min(rel_otu_neg_rel[rel_otu_neg_rel!=0]))
lowest


#pseudo = 4.702585e-06



kau_scores_probstel_rel<- igascores(posabunds = rel_otu_pos_rel, negabunds = 
rel_otu_neg_rel, pseudo =  4.702585e-06,  method = "kau")




######Steimle et al. Mild and Severe EAE mice groups###

#Mild EAE



otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_mice_ms_genus_mild_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_mice_ms_genus_mild_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_mice_ms_genus_mild_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_mice_ms_genus_mild_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


rownames(otu_pos_samples_df_inv) <- ss_mice_taxa$Genus
rownames(otu_neg_samples_df_inv) <- ss_mice_taxa$Genus


colnames(otu_pos_samples_df_inv) <- gsub("_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_mild <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_mild <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_mild[rel_otu_pos_mild!=0]),min(rel_otu_neg_mild[rel_otu_neg_mild!=0]))
lowest


#pseudo = 4.0803e-06



kau_scores_mild<- igascores(posabunds = rel_otu_pos_mild, negabunds = rel_otu_neg_mild, 
pseudo = 4.0803e-06,  method = "kau")


#Severe EAE


otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_mice_ms_genus_severe_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_mice_ms_genus_severe_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_mice_ms_genus_severe_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_mice_ms_genus_severe_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("_pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("_neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_severe <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_severe <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_severe[rel_otu_pos_severe!=0]),min(rel_otu_neg_severe[rel_otu_neg_severe!=0]))
lowest


#pseudo = 4.439472e-06



kau_scores_severe<- igascores(posabunds = rel_otu_pos_severe, negabunds = 
rel_otu_neg_severe, pseudo = 4.439472e-06,  method = "kau")






#Huus et al. Changes in IgA-targeted microbiota for reccurrent Clostridioides difficile 
infection
#CDI
#Pretransplant

otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_cdi_pre_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_cdi_pre_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_cdi_pre_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_cdi_pre_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_pre <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_pre <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_pre[rel_otu_pos_pre!=0]),min(rel_otu_neg_pre[rel_otu_neg_pre!=0]))
lowest


#pseudo = 1.324565e-05



kau_scores_cdi_pre<- igascores(posabunds = rel_otu_pos_pre, negabunds = rel_otu_neg_pre, 
pseudo = 1.324565e-05,  method = "kau")


kau_scores_cdi_pre_df<- rownames_to_column(kau_scores_cdi_pre, var = "Genus")



write_tsv(kau_scores_cdi_pre_df, 
"/Users/cimi_bioinformatics/Desktop/Bacteria_IgA_CIS/IC_score/kau_scores_cdi_pre_df.tsv")




##ProTransplant

otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_cdi_post_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_cdi_post_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_cdi_post_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_cdi_post_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_post <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_post <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_post[rel_otu_pos_post!=0]),min(rel_otu_neg_post[rel_otu_neg_post!=0]))
lowest


#pseudo =  1.419436e-05



kau_scores_cdi_post<- igascores(posabunds = rel_otu_pos_post, negabunds = rel_otu_neg_post, 
pseudo =  1.419436e-05,  method = "kau")



#Donor

otu_pos_samples_df <- as.data.frame(otu_table(phyloseq_cdi_donor_pos))
otu_pos_samples_df_inv <- as.data.frame(otu_pos_samples_df)


otu_neg_samples_df <- as.data.frame(otu_table(phyloseq_cdi_donor_neg))
otu_neg_samples_df_inv <- as.data.frame(otu_neg_samples_df)



rownames(otu_pos_samples_df_inv) <- NULL
rownames(otu_neg_samples_df_inv) <- NULL

tax_pos_samples <- as.data.frame(tax_table(phyloseq_cdi_donor_pos))
tax_neg_samples <- as.data.frame(tax_table(phyloseq_cdi_donor_neg))


rownames(otu_pos_samples_df_inv) <- make.names(tax_pos_samples$Genus, unique = TRUE)
rownames(otu_neg_samples_df_inv) <- make.names(tax_neg_samples$Genus, unique = TRUE)


colnames(otu_pos_samples_df_inv) <- gsub("pos", "", colnames(otu_pos_samples_df_inv))
colnames(otu_neg_samples_df_inv) <- gsub("neg", "", colnames(otu_neg_samples_df_inv))




otu_pos_samples_df_inv <- otu_pos_samples_df_inv[, order(colnames(otu_pos_samples_df_inv))]
otu_neg_samples_df_inv <- otu_neg_samples_df_inv[, order(colnames(otu_neg_samples_df_inv))]





rel_otu_pos_donor <- relabund(otu_pos_samples_df_inv)
rel_otu_neg_donor <- relabund(otu_neg_samples_df_inv)



lowest <- 
min(min(rel_otu_pos_donor[rel_otu_pos_donor!=0]),min(rel_otu_neg_donor[rel_otu_neg_donor!=0]))
lowest


#pseudo =  5.354466e-05



kau_scores_cdi_donor<- igascores(posabunds = rel_otu_pos_donor, negabunds = 
rel_otu_neg_donor, pseudo =  5.354466e-05,  method = "kau")




