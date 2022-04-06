library(reshape2)
library(ggplot2)
library(scales)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/SuperCorals/16S/")

map=read.table("Input_files/metadata.txt", header = T, sep = "\t", row.names = 1)
asv=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,1:16]
tax=read.table("outputs/SuperCorals_ASVs_noContanoOut.raw.txt", header = T, sep = "\t", row.names = 1)[,18:23]


## identify top 20 families 
names(asv)
fam.wid.agg=aggregate(asv, by = list(tax[, 4]), FUN =  sum)
topFamilies=fam.wid.agg[order(rowSums(fam.wid.agg[, 2:ncol(fam.wid.agg)]),decreasing = TRUE),][1:20,1]
fam.wid.agg$Group.1=ifelse(fam.wid.agg$Group.1 %in% topFamilies, as.character(fam.wid.agg$Group.1), "zOthers")
fa.gg2=aggregate(fam.wid.agg[, 2:ncol(fam.wid.agg)], by = list(fam.wid.agg[, 1]), FUN =  sum)
all.l=melt(fa.gg2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")

## Add sample information
all.l$Phenotype=map$Phenotype[match(all.l$Sample, rownames(map))]
all.l$Time=map$Time[match(all.l$Sample, rownames(map))]

final=all.l %>% group_by(Sample, Time, Phenotype, Family) %>% dplyr::summarise(Abundance=sum(Abundance))

## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
#P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")

pdf("./outputs/replicated_barplots_families.pdf",  width = 6, height =5, pointsize = 12) 
ggplot() + geom_bar(aes(y = Abundance, x = Sample, fill = Family), 
                   data = final, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + 
  scale_fill_manual(values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), 
         legend.position = 'right', axis.text.x=element_text(angle=90,hjust=1)) + 
  guides(fill=guide_legend(ncol=1))  +
  facet_wrap(Phenotype~Time,nrow =2,scales = "free_x") 
dev.off() 


## identify top 20 genera
names(asv)
fam.wid.agg=aggregate(asv, by = list(tax[, 5]), FUN =  sum)
topFamilies=fam.wid.agg[order(rowSums(fam.wid.agg[, 2:ncol(fam.wid.agg)]),decreasing = TRUE),][1:20,1]
fam.wid.agg$Group.1=ifelse(fam.wid.agg$Group.1 %in% topFamilies, as.character(fam.wid.agg$Group.1), "zOthers")
fa.gg2=aggregate(fam.wid.agg[, 2:ncol(fam.wid.agg)], by = list(fam.wid.agg[, 1]), FUN =  sum)
all.l=melt(fa.gg2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Genus","Sample","Abundance")

## Add sample information
all.l$Phenotype=map$Phenotype[match(all.l$Sample, rownames(map))]
all.l$Time=map$Time[match(all.l$Sample, rownames(map))]

final=all.l %>% group_by(Sample, Time, Phenotype, Genus) %>% dplyr::summarise(Abundance=sum(Abundance))

## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
#P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")

pdf("./outputs/replicated_barplots_genus.pdf",  width = 6, height =5, pointsize = 12) 
ggplot() + geom_bar(aes(y = Abundance, x = Sample, fill = Genus), 
                    data = final, stat="identity", position = "fill")  +
  labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + 
  scale_fill_manual(values=P21)  +   theme_classic() + 
  theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), 
         legend.position = 'right', axis.text.x=element_text(angle=90,hjust=1)) + 
  guides(fill=guide_legend(ncol=1))  +
  facet_wrap(Phenotype~Time,nrow =2,scales = "free_x") 
dev.off() 
