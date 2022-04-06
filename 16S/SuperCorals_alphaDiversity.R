library(vegan)
library(ggplot2)
library(GUniFrac)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/SuperCorals/16S/")

map=read.table("Input_files/metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,1:16]
tax=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,18:23]

###rarefying
cnts=t(asv)
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, 8582)$otu.tab.rff

############################################################
##################### Alpha-diversity ######################
############################################################

alpha=as.data.frame(t(estimateR(asv.rar)))
alpha$Shannon=diversity(asv.rar, index = "shannon")
alpha$Time=map$Time[match(rownames(alpha), rownames(map))]
alpha$Phenotype=map$Phenotype[match(rownames(alpha), rownames(map))]

##################################################
##################### Stats ######################
##################################################
#Shannon
shapiro.test(alpha$Shannon) # p-value > 0.05 implying we can assume normality.
anova_shan=summary(aov(alpha$index ~ alpha$Phenotype * alpha$Time))
anova_shan_df=as.data.frame(anova[[1]])
write.table(anova_shan_df, "outputs/anova_shannon.txt", sep = "\t", row.names = T, quote = F)
TukeyHSD(aov(alpha$Shannon ~ alpha$Phenotype * alpha$Time))


phen_t1=subset(alpha, Time == "F1" )
pairwise.t.test(phen_t1$Shannon, phen_t1$Phenotype,  p.adj = "fdr", digits = 5)[[3]]
phen_t2=subset(alpha, Time == "F2" )
pairwise.t.test(phen_t2$Shannon, phen_t2$Phenotype,  p.adj = "fdr", digits = 5)[[3]]

sen=subset(alpha, Phenotype == "Sensitive" )
pairwise.t.test(sen$Shannon, sen$Time,  p.adj = "fdr", digits = 5)[[3]]
res=subset(alpha, Phenotype == "Resistant" )
pairwise.t.test(res$Shannon, res$Time,  p.adj = "fdr", digits = 5)[[3]]

## Chao1
shapiro.test(alpha$S.chao1) # p-value > 0.05 implying we can assume normality.
anova_cha=summary(aov(alpha$S.chao1 ~ alpha$Phenotype * alpha$Time))
anova_cha_df=as.data.frame(anova_cha[[1]])
write.table(anova_cha_df, "outputs/anova_chao.txt", sep = "\t", row.names = T, quote = F)
TukeyHSD(aov(alpha$S.chao1 ~ alpha$Phenotype * alpha$Time))

phen_t1=subset(alpha, Time == "F1" )
pairwise.t.test(phen_t1$S.chao1, phen_t1$Phenotype,  p.adj = "fdr", digits = 5)[[3]]
phen_t2=subset(alpha, Time == "F2" )
pairwise.t.test(phen_t2$S.chao1, phen_t2$Phenotype,  p.adj = "fdr", digits = 5)[[3]]

sen=subset(alpha, Phenotype == "Sensitive" )
pairwise.t.test(sen$Shannon, sen$Time,  p.adj = "fdr", digits = 5)[[3]]
res=subset(alpha, Phenotype == "Resistant" )
pairwise.t.test(res$Shannon, res$Time,  p.adj = "fdr", digits = 5)[[3]]


# boxplots
shan=ggplot(alpha, aes(x=Phenotype, y=Shannon, fill=Phenotype)) + 
  stat_boxplot(geom = "errorbar")  + 
  geom_boxplot(alpha = 1) +  
  scale_fill_manual(values=c("#00607A", "#A44200"))  + 
  facet_grid(~Time) + 
  theme_bw() + labs( y= "Shannon diversity", x="") 


cha1=ggplot(alpha, aes(x=Phenotype, y=S.chao1, fill=Phenotype)) + 
  stat_boxplot(geom = "errorbar")  + 
  geom_boxplot(alpha = 1) +  
  scale_fill_manual(values=c("#00607A", "#A44200"))  + 
  facet_grid(~Time) + 
  theme_bw() + labs( y= "Chao1 estimated richness", x="") 

pdf("./outputs/SuperCorals_AlphaDiversiy.pdf", width=5,height=5, pointsize = 12)
shan/cha1
dev.off()
