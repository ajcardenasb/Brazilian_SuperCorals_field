library(ANCOMBC)
library(phyloseq)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/SuperCorals/16S/")

map=read.table("Input_files/metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,1:16]
tax=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,18:23]
otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)
phy.gen=aggregate_taxa(phy, "Genus")

##################################
## between phenotypes at T1|T2 ##
#################################
names(map)
phy_phen_t1=subset_samples(phy.gen, Time == "F1")
res1=ancombc(phyloseq=phy_phen_t1,formula="Phenotype",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Phenotype",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res1_df=data.frame(res1[["res"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "Resistant", "Sensitive")
res1_sig$Genus=tax$Genus[match(rownames(res1_sig),tax$Genus)]
res1_sig$Family=tax$Family[match(rownames(res1_sig),tax$Genus)]
res1_sig$Comparison="RxS F1"
message("Number of total DA genera: ", nrow(res1_sig), "\nNumber of DA genera enriched in Resistant: ", nrow(subset(res1_sig, Diff_more_abundant == "Resistant" )),  "\nNumber of DA genera enriched in Sensitive: ", nrow(subset(res1_sig, Diff_more_abundant == "Sensitive")))

phy_phen_t2=subset_samples(phy.gen, Time == "F2")
res2=ancombc(phyloseq=phy_phen_t2,formula="Phenotype",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Phenotype",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res2_df=data.frame(res2[["res"]])
colnames(res2_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "Resistant", "Sensitive")
res2_sig$Genus=tax$Genus[match(rownames(res2_sig),tax$Genus)]
res2_sig$Family=tax$Family[match(rownames(res2_sig),tax$Genus)]
res2_sig$Comparison="RxS F2"
message("Number of total DA genera: ", nrow(res2_sig), "\nNumber of DA genera enriched in Resistant: ", nrow(subset(res2_sig, Diff_more_abundant == "Resistant" )),  "\nNumber of DA genera enriched in Sensitive: ", nrow(subset(res2_sig, Diff_more_abundant == "Sensitive")))

##########################
## between timepoints  ##
#########################

sen=subset_samples(phy.gen, Phenotype == "Sensitive")
res3=ancombc(phyloseq=sen,formula="Time",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Time",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res3_df=data.frame(res3[["res"]])
colnames(res3_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "F1", "F2")
res3_sig$Genus=tax$Genus[match(rownames(res3_sig),tax$Genus)]
res3_sig$Family=tax$Family[match(rownames(res3_sig),tax$Genus)]
res3_sig$Comparison="Sen F1xF2"
message("Number of total DA genera: ", nrow(res3_sig), "\nNumber of DA genera enriched in F1: ", nrow(subset(res3_sig, Diff_more_abundant == "F1" )),  "\nNumber of DA genera enriched in F2: ", nrow(subset(res3_sig, Diff_more_abundant == "F2")))

res=subset_samples(phy.gen, Phenotype == "Resistant")
res4=ancombc(phyloseq=res,formula="Time",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Time",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = T,alpha = 0.05,global = TRUE)
res4_df=data.frame(res4[["res"]])
colnames(res4_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res4_sig=subset(res4_df, Diff_abundant == "TRUE")
res4_sig$Diff_more_abundant=ifelse( res4_sig$W < 0 , "F1", "F2")
res4_sig$Genus=tax$Genus[match(rownames(res4_sig),tax$Genus)]
res4_sig$Family=tax$Family[match(rownames(res4_sig),tax$Genus)]
res4_sig$Comparison="Res F1xF2"
message("Number of total DA genera: ", nrow(res4_sig), "\nNumber of DA genera enriched in F1: ", nrow(subset(res4_sig, Diff_more_abundant == "F1" )),  "\nNumber of DA genera enriched in F2: ", nrow(subset(res4_sig, Diff_more_abundant == "F2")))

all=rbind(res1_sig, res2_sig, res3_sig, res4_sig)
