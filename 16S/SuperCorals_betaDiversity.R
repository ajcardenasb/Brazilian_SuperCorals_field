library(phyloseq)
library(ggplot2)
library(vegan)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/SuperCorals/16S/")

map=read.table("Input_files/metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,1:16]
tax=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,18:23]
otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t=tax_table(as.matrix(tax))
phy= phyloseq(otu.t,  sam.t , tax.t)

P2=c("#00607A", "#A44200")

phy.t=microbiome::transform(phy, transform = "clr", target = "OTU", shift = 0, scale = 1)

ord = ordinate(phy.t, method = "RDA", distance = "euclidean")
pdf("./outputs/16S_ordination.pdf", width=4.5,height=3, pointsize = 12)
plot_ordination(phy.t,ord, color = "Phenotype", shape = "Time") + 
  geom_point(size = 4, alpha = 1)  + 
  scale_colour_manual(values=P2) + ggtitle("") +  
  theme_bw() + theme( legend.position = 'bottom')
dev.off()


#####################################################
######## Stats on community composition ############
####################################################
library(vegan)
library(pairwiseAdonis)
#otu.n=as.data.frame(t(sweep(otu,2,colSums(otu),"/")))
otu.n=data.frame(t(otu_table(phy.t)))
otu.n$Time=map$Time[match(rownames(otu.n), rownames(map))]
otu.n$Phenotype=map$Phenotype[match(rownames(otu.n), rownames(map))]
otu.n$group=paste(otu.n$Phenotype,otu.n$Time)

#betadisper
distance=vegdist(otu.n[,1:3520], method = "euclidean")

beta_phenotype=betadisper(distance, otu.n$Phenotype)
permutest(beta_phenotype, permutations = 99, pairwise = TRUE) 

beta_time=betadisper(distance, otu.n$Time)
permutest(beta_time, permutations = 99, pairwise = TRUE) 

beta=betadisper(distance, otu.n$group)
permutest(beta, permutations = 99, pairwise = TRUE) 

##overal model
adonis=adonis(otu.n[,1:3520]~ otu.n$Phenotype * otu.n$Time, method = "euclidean" )
adonis_df=as.data.frame(adonis[["aov.tab"]])
write.table(adonis_df, "outputs/overall_adonis.txt", sep = "\t", row.names = T, quote = F)

#all comparisons
all_pairWadonis_df=pairwise.adonis(otu.n[,1:3520], otu.n$group,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) 
write.table(species_pairWadonis_df, "outputs/pairwiseAdonis_field.txt", sep = "\t", row.names = F, quote = F)
