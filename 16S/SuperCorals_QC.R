setwd("~/Documents/Bioinformatics_scripts/R_scripts/SuperCorals/16S/")
library(stringr)


##################################################################
#####Identifying and removing contaminant ASVs normalized data#####
###################################################################

asv=read.table("Input_files/SuperCorals_ASV_table_pooled.txt", header = TRUE, row.names = 1, sep = " ")[,1:17]
tax=read.table("Input_files/SuperCorals_ASV_table_pooled.txt", header = TRUE, row.names = 1, sep = " ")[,19:25]
#map=read.table("Input_files/metadata.txt", header = T, row.names = 1, sep = "\t")

# all samples have more than 10000 reads
hist(colSums(asv),  breaks = 50, labels = F)
axis(side=1, at=seq(0,150000, 10000))

asv.o=asv[, colSums(asv) > 10000]
message(ncol(asv.o)," samples with > 1000 reads were retained out of ", ncol(asv), " total samples")

#Identify and removing contaminant ASVs raw data
asv.r=as.data.frame(sweep(asv.o,2,colSums(asv.o),`/`))
asv.r$Sum=rowSums(asv.r[,2:ncol(asv.r)])
names(asv.r)
asv.r$contaFactor=(asv.r$BlankSample/asv.r$Sum)*100
rownames(asv.r)=rownames(asv)
Conta=subset(asv.r, asv.r$contaFactor > 10)
Conta$Family=tax$Family[match(rownames(Conta), rownames(asv))]
message("Number of total ASVs: ", nrow(asv))
message("Number of identified contaminant ASVs removed from the analysis: ", length(rownames(Conta)), "\n", Conta$Family[1],"\n", Conta$Family[2],"\n", Conta$Family[3],"\n", Conta$Family[4],"\n", Conta$Family[5])

#remove any chloroplast or mitochobdria
unwant_tax=tax %>% filter_all(any_vars(str_detect(., 'Mitochondria|Chloroplast')))
colnames(asv.o)
asv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(Conta) & !rownames(asv.o) %in% rownames(unwant_tax))[,-c(1)]
colnames(asv.noConta)
dim(asv.noConta)

#remove ASVs with only zeros
asv.noRare=asv.noConta[rowSums(asv.noConta[])>0,]

# Export normalized and raw ASV tables
asv.noConta.f=merge(asv.noRare, tax, by="row.names")

#rename sample name
write.table(asv.noConta.f, "./outputs/SuperCorals_ASVs_noContanoOut.raw.txt",  quote = FALSE, row.names=F, sep = "\t") 
message("Number of ASVs used in the analysis: ", length(rownames(asv.noConta.f)))

