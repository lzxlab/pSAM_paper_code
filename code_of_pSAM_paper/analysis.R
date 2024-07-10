setwd("/home/yukai/projects/CBAmodel/subCellLoc/bin")
library(data.table)
library(clusterProfiler)
library(dplyr)
library(msigdbr)
library(ggplot2)
library(scales)
library(ggrepel)
library(magrittr)
library(ggpubr)
library(stringr)
library(survival)
library(survminer)
library(matrixStats)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
library(limma)
library(pROC)
library(tidyr)
#################################functions#################################
survplot <- function(dat = mydata, type = "OS", fit = fit, pval = pval){
  p <- ggsurvplot(fit,
                  linetype = 1,
                  #censor.shape=45,
                  data = dat,
                  size = 1, # change line size
                  palette = c("#02BFC4", "#F8766D"),# custom color palettes
                  #conf.int = TRUE, # Add confidence interval
                  pval = paste('p = ', round(pval, 3)), # Add p-value
                  risk.table = T, # Add risk table
                  #tables.theme = theme_survminer(font.main = 10),
                  #risk.table.col = "strata",# Risk table color by groups
                  legend = "right",
                  #legend.labs = c("G1 (n = 7)", "G2 (n = 26)", "G3 (n = 24)"), # Change legend labels
                  risk.table.height = 0.25, # Useful to change when you have multiple groups
                  ggtheme = theme_classic2(), # Change ggplot2 theme
                  xlab = "Time (days)",
                  ylab = paste0("Probability of ", type))
  return(p)
}
str2lst <- function(ling, fullength){
  tmp_list = c()
  tmp_str = strsplit(ling, ";")
  for (i in tmp_str[[1]]) {
    tmp_i = strsplit(i, "[.][.][.]")
    tmp_list = c(tmp_list, c(max(as.integer(tmp_i[[1]][1]), 1):min(as.integer(tmp_i[[1]][2]), fullength)))
  }
  return(unique(tmp_list))
}
str2lst2 <- function(ling, fullength){
  tmp_list = c()
  tmp_str = strsplit(ling, ";")
  for (i in tmp_str[[1]]) {
    tmp_i = strsplit(i, "[.][.][.]")
    tmp_list = c(tmp_list, c(max(as.integer(tmp_i[[1]][1]) - 16, 1):min(as.integer(tmp_i[[1]][2]) + 16, fullength)))
  }
  return(unique(tmp_list))
}
readgmt <- function(gmtfile){
  kegggmt <- strsplit(gmtfile, "\t")
  names(kegggmt) <- vapply(kegggmt, function(y) y[1], character(1))
  kegggmt <- lapply(kegggmt, "[", -c(1:2))
  
  pathways <- c()
  genes <- c()
  for (i in 1:length(kegggmt)) {
    if (length(unlist(kegggmt[i])) > 0){
      for (j in 1:length(unlist(kegggmt[i]))){
        pathways <- c(pathways, names(kegggmt[i]))
        genes <- c(genes, unlist(kegggmt[i])[j])
      }
    }
    print(paste(i, " of ", length(kegggmt), " has finished!"))
  }
  
  keggpath2gene <- data.frame(gs_name = pathways,
                              genes = as.character(genes))
  return(keggpath2gene)
}
two_matrix_cor <- function(x, y){
  if (!is.matrix(x)) x <- matrix(x, ncol=1L)
  if (!is.matrix(y)) y <- matrix(y, ncol=1L)
  ncx <- ncol(x)
  ncy <- ncol(y)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  p <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    print(paste('mat1 in line ', i, ' of all ', ncx, ' samples!!!!', sep = ' '))
    for (j in seq_len(ncy)) {
      #print(paste('mat1 in line ', i, ' of all ', j, ' samples!!!!', sep = ' '))
      x2 <- x[,i]
      y2 <- y[,j]
      ok <- complete.cases(x2, y2)
      x2 <- x2[ok]
      y2 <- y2[ok]
      r[i, j] <- if(length(ok[ok]) > 2) as.numeric(cor.test(x2, y2)$estimate) else NA
      p[i, j] <- if(length(ok[ok]) > 2) as.numeric(cor.test(x2, y2)$p.value) else NA
    }
  }
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(y)
  rownames(p) <- colnames(x)
  colnames(p) <- colnames(y)
  cor_r <- melt(r)
  cor_p <- melt(p)
  cor_res <- cbind(cor_r, cor_p[, 3])
  names(cor_res) <- c("ID1", "ID2", "cor", "pva")
  return(cor_res)
}
surv_bestcut <- function(mtr, gene, cli_info, num = 20){
  tmp_mtr <- mtr[which(rownames(mtr) == gene), ]
  if (length(tmp_mtr[is.na(tmp_mtr)]) == 0) {
    tmp_mtr <- tmp_mtr
  }else{
    tmp_mtr <- tmp_mtr[-which(is.na(tmp_mtr))]
  }
  common_samples <- intersect(names(tmp_mtr), rownames(cli_info))
  cluster_surv <- cli_info[common_samples, ]
  tmp_mtr <- as.data.frame(t(tmp_mtr))[common_samples, ]
  sevalue <- as.numeric(tmp_mtr)
  values <- c()
  hr_os <- c()
  pva_os <- c()
  hr_rfs <- c()
  pva_rfs <- c()
  n_high <- c()
  n_low <- c()
  for (i in c(round(length(sevalue)/4):round(length(sevalue)/4*3))) {
    cluster_surv$Type = ifelse(sevalue > sort(sevalue)[i], "0.High", "1.Low")
    values <- c(values, sort(sevalue)[i])
    tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
    hr_os <- c(hr_os, tmp$conf.int[[1]])
    pva_os <- c(pva_os, tmp$logtest[[3]])
    ##PFS
    #tmp <- summary(coxph((Surv(DFI.time, DFI)) ~ Type, data = cluster_surv))
    #hr_rfs <- c(hr_rfs, tmp$conf.int[[1]])
    #pva_rfs <- c(pva_rfs, tmp$logtest[[3]])
    n_high <- c(n_high, (length(sevalue)-i))
    n_low <- c(n_low, i)
  }
  res <- data.frame(ID = gene,
                    cutoff = values,
                    HR_OS = hr_os,
                    Pvalue_OS = pva_os,
                    n_high = n_high,
                    n_low = n_low)
  
}

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(face = "bold",
                            size = 18, color = 'black', hjust = 0.5),
    axis.title = element_text(face = "bold",
                              size = 15, color = 'black')
  )
#################################architecture of paper#################################
##mutation, subloc, DNN, tumor, pancancer
#################################effect of mutation to subloc#################################
#mito:      TargetP2.0, Mito453
#membrane:    
#nucleus:   NLS (SeqNLS, UniProt, NLSdb), NES (ValidNESs, NESbase, UniProt, NLSdb) 
#cyto:   
#secreted:  Signal Peptide Website, TargetP2.0
##NLS, NES, SP, MT
########mutation enriched in signals for each type########
###NLS NES SP MT
#NLS
sigregions <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/inuloc.point.mut.NLSregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#NES
sigregions <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/inuloc.point.mut.NESregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#SP
sigregions <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/inuloc.point.mut.SPregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)

#MT
sigregions <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/inuloc.point.mut.MTregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)


##plot
res <- fread("../8.results/inuloc.point.mut.NLSregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot1 <- melt(res, id.vars = c("ID", "Tumor"))
p1 <- ggplot(res4plot1, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (NLS)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/inuloc.point.mut.NESregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot2 <- melt(res, id.vars = c("ID", "Tumor"))
p2 <- ggplot(res4plot2, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (NES)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/inuloc.point.mut.SPregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot3 <- melt(res, id.vars = c("ID", "Tumor"))
p3 <- ggplot(res4plot3, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (SP)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/inuloc.point.mut.MTregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot4 <- melt(res, id.vars = c("ID", "Tumor"))
p4 <- ggplot(res4plot4, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (MT)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("../8.results/inuloc.point.mut.distibute.conbined.boxplot.pdf", width = 10, height = 16)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()

########mutation enriched in signals for combined########
sigregions1 <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions2 <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions3 <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions4 <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions <- rbind(sigregions1, sigregions2, sigregions3, sigregions4)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- c("Pancancer")
res <- data.frame()
for (cancer in tumorlist) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/inuloc.point.mut.Mergedregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)

res <- fread("../8.results/inuloc.point.mut.Mergedregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot4 <- melt(res, id.vars = c("ID", "Tumor"))
p1 <- ggplot(res4plot4, aes(x = variable, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (MT)")+
  theme_classic2()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("../8.results/inuloc.point.mut.distibute.Merged.boxplot.pdf", width = 4, height = 4)
print(p1)
dev.off()

########annote mutation file########
mutfile <- fread("/home/yukai/data/TCGA/mutect_tcga/TCGA-BLCA.mutect2_snv.tsv",
                 header = T, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf", header = F, stringsAsFactors = F, data.table = F)
names(singmut) <- c(names(mutfile)[2:(ncol(mutfile) -1)], "UniID")
annovarscore <- fread("/home/yukai/data/ANNOVAR/hg38_dbnsfp30a.txt",
                      header = T, stringsAsFactors = F, data.table = F,
                      select = c("#Chr", "Start", "Ref", "Alt", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score",
                                 "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
                                 "LRT_score", "LRT_pred", "FATHMM_score", "FATHMM_pred"))
names(annovarscore)[1] <- "Chr"
annovarscore$ID <- paste(annovarscore$Chr, annovarscore$Start, annovarscore$Ref, annovarscore$Alt, sep="_")
singmut$ID <- paste(gsub("chr", "", singmut$chrom), singmut$start, singmut$ref, singmut$alt, sep="_")
singmut <- merge(singmut, annovarscore, by.x = "ID", by.y = "ID")
write.table(singmut, "../4.panMut/merged.mut.maf.annote",
            row.names = F, col.names = T, sep = "\t", quote = F)

########mutation deleterous for each type########
sigregions <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf.annote", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  pro <- as.character(res[i, ]$UniID)
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}
res$Type <- genes
write.table(res, "../8.results/inuloc.point.mut.NLSdeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf.annote", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  pro <- as.character(res[i, ]$UniID)
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}
res$Type <- genes
write.table(res, "../8.results/inuloc.point.mut.NESdeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf.annote", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  pro <- as.character(res[i, ]$UniID)
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}
res$Type <- genes
write.table(res, "../8.results/inuloc.point.mut.SPdeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf.annote", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  pro <- as.character(res[i, ]$UniID)
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}
res$Type <- genes
write.table(res, "../8.results/inuloc.point.mut.MTdeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)

##plot
res <- fread("../8.results/inuloc.point.mut.NLSdeleterous.txt", header = T, stringsAsFactors = F, data.table = F)
res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)
tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("0.NLS", "1.Out", "0.NLS", "1.Out"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)
sigif <- prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p1 <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab(sigif$p.value)+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()

res <- fread("../8.results/inuloc.point.mut.NESdeleterous.txt", header = T, stringsAsFactors = F, data.table = F)
res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)
tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("1.Out", "0.NES", "1.Out", "0.NES"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)
sigif <- prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p2 <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab(sigif$p.value)+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()

res <- fread("../8.results/inuloc.point.mut.SPdeleterous.txt", header = T, stringsAsFactors = F, data.table = F)
res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)
tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("1.Out", "0.SP", "1.Out", "0.SP"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)
sigif <- prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p3 <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab(sigif$p.value)+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()

res <- fread("../8.results/inuloc.point.mut.MTdeleterous.txt", header = T, stringsAsFactors = F, data.table = F)
res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)
tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("1.Out", "0.MT", "1.Out", "0.MT"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)
sigif <- prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p4 <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab(sigif$p.value)+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()

pdf("../8.results/inuloc.point.mut.deleterous.eachtype.pdf", width = 8, height =5)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()


########mutation deleterous for combined########
sigregions1 <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions2 <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions3 <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions4 <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions <- rbind(sigregions1, sigregions2, sigregions3, sigregions4)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
singmut <- fread("../4.panMut/merged.mut.maf.annote", header = T, stringsAsFactors = F, data.table = F)
singmut$Merge1 <- ifelse(singmut$SIFT_pred == "D", 1, 0)
singmut$Merge2 <- ifelse(singmut$Polyphen2_HDIV_pred == "D", 1, 0)
singmut$Merge3 <- ifelse(singmut$Polyphen2_HVAR_pred == "D", 1, 0)
singmut$Merge4 <- ifelse(singmut$LRT_pred == "D", 1, 0)
singmut$Merge5 <- ifelse(singmut$FATHMM_pred == "D", 1, 0)
singmut$Level <- singmut$Merge1 + singmut$Merge2 + singmut$Merge3 + singmut$Merge4 + singmut$Merge5
res <- singmut[, c(2, 8, 11, ncol(singmut))]
res <- res[res$UniID %in% sigregions$V1, ]
genes <- c()
num = 1
for (i in c(1:nrow(res))) {
  tmpmut <- res[i, ]
  pro <- as.character(res[i, ]$UniID)
  mutpos <- as.numeric(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")[2])
  proseq <- c(1:nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
  if (mutpos %in% nlsregion) {
    genes <- c(genes, "NLS")
  }else{
    genes <- c(genes, "Out")
  }
  print(paste(num, " of ", nrow(res), " has been finished!!!"))
  num = num + 1
}
res$Type <- genes
write.table(res, "../8.results/inuloc.point.mut.Mergeddeleterous.txt", col.names = T, row.names = F, sep = '\t', quote = F)

res <- fread("../8.results/inuloc.point.mut.Mergeddeleterous.txt", header = T, stringsAsFactors = F, data.table = F)
res4plot <- as.data.frame.array(table(res[, c(4,5)]))
res4plot$NLSRate <- res4plot$NLS/sum(res4plot$NLS)
res4plot$OutRate <- res4plot$Out/sum(res4plot$Out)
res4plot$Type <- rownames(res4plot)
#wilcox.test(res[res$Type == "NLS", ]$Level, res[res$Type != "NLS", ]$Level)
tmpplot <- data.frame(ID = c("0-2", "0-2", "3-5", "3-5"),
                      Type = c("1.Out", "0.TPs", "1.Out", "0.TPs"),
                      Rates = c(sum(res4plot[res4plot$Type %in% c(0,1,2), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(0,1,2), ]$OutRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$NLSRate),
                                sum(res4plot[res4plot$Type %in% c(3,4,5), ]$OutRate)),
                      stringsAsFactors = F)
sigif <- prop.test(c(sum(res4plot$NLS[c(1,2,3)]), sum(res4plot$Out[c(1,2,3)])), c(sum(res4plot$NLS), sum(res4plot$Out)))
p4 <- ggplot(tmpplot, aes(x = ID, y = Rates, fill = Type))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab(sigif$p.value)+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Relative Rates")+
  theme_classic2()

pdf("../8.results/inuloc.point.mut.deleterous.Merged.pdf", width = 4, height =2.5)
print(p4)
dev.off()

########classification of subloc########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]

res <- as.data.frame(table(subloc$V3))
res$Group <- paste(res$Var1, "( N = ", res$Freq, " )", sep = "")
clinical_table <- merge(subloc, res, by.x = "V3", by.y = "Var1")
colorinfo <- brewer.pal(5, "Set1")
##Age
p1 <- ggplot(data=clinical_table, mapping=aes(x="Group",fill=Group))+
  geom_bar(stat="count",width=0.5,position='stack',size=5)+
  coord_polar("y", start=0)+
  scale_fill_manual(values=colorinfo)+
  ggtitle('Class')+
  blank_theme +
  geom_text(stat="count",aes(label = scales::percent(..count../length(rownames(clinical_table)))), size=4, position=position_stack(vjust = 0.5))

p1
ggsave("../8.results/inuloc.sublocs.distribut.pdf", p1, width = 6, height = 6)

########single and multi loc distribution########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
colors <- c("#984EA3", "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- unique(subloc[, c(1, 3)])

res <- as.data.frame(table(subloc$V1))
res4plot <- as.data.frame(table(res$Freq))
res4plot$Var1 <- factor(res4plot$Var1, c(1,2,3,4,5))
res4plot$Freq <- log10(res4plot$Freq)
p4 <- ggplot(res4plot, aes(x = Var1, y = Freq, fill = Var1))+
  geom_bar(stat="identity", position=position_dodge())+
  #stat_compare_means(label = "p.signif")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("log10 Protein Counts")+
  theme_classic2()

pdf("../8.results/inuloc.sublocs.distribut.multiloc.pdf", width = 4, height =2.5)
print(p4)
dev.off()

########multiloc pros oncoscore########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- unique(subloc[, c(1, 3)])

res <- as.data.frame(table(subloc$V1))
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
res <- merge(res, hsdataset, by.x = "Var1", by.y = "V1")

e3s <- res[res$Freq > 2, ]$V2
allfiles <- fread("/home/yukai/softwares/2020plus-master/tumor_sample_stat.txt", header = F, stringsAsFactors = F, data.table = F)
aikeids <- allfiles$V1

finalres <- data.frame()
for (ids in aikeids) {
  scores <- fread(paste("/home/yukai/softwares/2020plus-master/output_res/", ids, "/output/results/r_random_forest_prediction.txt", sep = ""), header = T, stringsAsFactors = 1, data.table = F)
  scores <- scores[, c(1, 26, 27, 37, 38)]
  scores$Class <- ids
  finalres <- rbind(finalres, scores)
}

##uniq and densityplot
res4plot <- finalres
res4plot <- res4plot[res4plot$gene %in% e3s, ]
res4plot1 <- res4plot[order(res4plot$oncogene, decreasing = T), ]
res4plot1 <- res4plot1 %>% distinct(gene, .keep_all = TRUE)
res4plot2 <- res4plot[order(res4plot$tsg, decreasing = T), ]
res4plot2 <- res4plot2 %>% distinct(gene, .keep_all = TRUE)
res4plot3 <- res4plot[order(res4plot$`driver p-value`), ]
res4plot3 <- res4plot3 %>% distinct(gene, .keep_all = TRUE)
res4plot1 <- res4plot1[order(res4plot1$gene), ]
res4plot2 <- res4plot2[order(res4plot2$gene), ]
res4plot3 <- res4plot3[order(res4plot3$gene), ]
res4plot <- data.frame(gene = res4plot1$gene,
                       oncogene = res4plot1$`oncogene score`,
                       tsg = res4plot2$`tsg score`,
                       pvalue = res4plot3$`driver p-value`)

res4plot$group <- ifelse(res4plot$pvalue < 0.05, "Significant", "None")
res4plot4sig <- res4plot[res4plot$pvalue < 0.005, ]
p1 <- ggplot(res4plot, aes(x = tsg, y = oncogene))+
  geom_point(aes(x = tsg, y = oncogene, color = group), size = 0.6)+
  stat_density_2d(colour="black", lwd = 0.4)+
  scale_color_manual(values = c("gray", "red"))+
  xlab("TSG score")+
  ylab("Oncogene score")+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  geom_text_repel(data = res4plot4sig, aes(label = gene), color = "red", size = 2)+
  theme_bw()
p1
ggsave("../8.results/inuloc.multiloc.genes.oncoscores.pdf", p1, width = 4, height = 2.8)

########multiloc pros CERES########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- unique(subloc[, c(1, 3)])

res <- as.data.frame(table(subloc$V1))
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
res <- merge(res, hsdataset, by.x = "Var1", by.y = "V1")

e3s <- res[res$Freq > 2, ]$V2
allfiles <- fread("/home/yukai/softwares/2020plus-master/tumor_sample_stat.txt", header = F, stringsAsFactors = F, data.table = F)
aikeids <- allfiles$V1
finalres <- data.frame()
for (ids in aikeids) {
  scores <- fread(paste("/home/yukai/softwares/2020plus-master/output_res/", ids, "/output/results/r_random_forest_prediction.txt", sep = ""), header = T, stringsAsFactors = 1, data.table = F)
  scores <- scores[, c(1, 26, 27, 37, 38)]
  scores$Class <- ids
  finalres <- rbind(finalres, scores)
}

##uniq and densityplot
res4plot <- finalres
res4plot <- res4plot[res4plot$gene %in% e3s, ]
res4plot1 <- res4plot[order(res4plot$oncogene, decreasing = T), ]
res4plot1 <- res4plot1 %>% distinct(gene, .keep_all = TRUE)
res4plot2 <- res4plot[order(res4plot$tsg, decreasing = T), ]
res4plot2 <- res4plot2 %>% distinct(gene, .keep_all = TRUE)
res4plot3 <- res4plot[order(res4plot$`driver p-value`), ]
res4plot3 <- res4plot3 %>% distinct(gene, .keep_all = TRUE)
res4plot1 <- res4plot1[order(res4plot1$gene), ]
res4plot2 <- res4plot2[order(res4plot2$gene), ]
res4plot3 <- res4plot3[order(res4plot3$gene), ]
res4plot <- data.frame(gene = res4plot1$gene,
                       oncogene = res4plot1$`oncogene score`,
                       tsg = res4plot2$`tsg score`,
                       pvalue = res4plot3$`driver p-value`)

res4plot$group <- ifelse(res4plot$pvalue < 0.01, "Significant", "None")
siggenes <- unique(res4plot[res4plot$group == "Significant", ]$gene)

ceres <- fread("/home/yukai/data/Datasets/DeMap/CCLE_D2_combined_gene_dep_scores.csv", header = T, stringsAsFactors = F, data.table = F)
ceres$V1 <- gsub(" .*", "", ceres$V1)
rownames(ceres) <- ceres$V1

alltiss <- as.data.frame(table(sub("[^_]*_", "", names(ceres))))
alltiss <- alltiss[alltiss$Freq > 20, ]

e3ceres <- ceres[intersect(rownames(ceres), e3s), ]
e3sval <- data.frame(ID = rownames(e3ceres),
                     Median = rowMedians(as.matrix(e3ceres[, -1]), na.rm = T))
e3sval <- e3sval[order(e3sval$Median), ]

tmpceres <- ceres[, -1]
ceres_md <- data.frame(ID = rownames(tmpceres),
                       Score = rowMeans(tmpceres, na.rm = T),
                       Class = "All")
wca <- ceres_md
wca$Class <- ifelse(wca$ID %in% e3sval$ID[1:30], "0.Sigpros", "1.Others")
wca <- wca[is.nan(wca$Score) == FALSE, ]
library(plyr)
mu <- ddply(wca, "Class", summarise, grp.mean=mean(Score))
p<-ggplot(wca, aes(x=Score, fill=Class)) +
  geom_density(alpha=0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Class),
             linetype="dashed")+
  theme_classic2()
print(p)
ggsave("../8.results/inuloc.multiloc.genes.ceres.pdf", p, width = 6, height = 5)


########multiloc pros enrichment########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- unique(subloc[, c(1, 3)])

res <- as.data.frame(table(subloc$V1))
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
res <- merge(res, hsdataset, by.x = "Var1", by.y = "V1")

e3s <- res[res$Freq > 2, ]$V2

kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

res_hall <- enricher(e3s, pvalueCutoff = 1,
                     TERM2GENE = kegg, qvalueCutoff = 1, minGSSize = 1)
res2subtop15 <- res_hall@result
res2subtop15$link <- -log10(res2subtop15$p.adjust)
res2subtop15 <- res2subtop15[res2subtop15$p.adjust < 0.00001, ]
res2subtop15$ID <- factor(res2subtop15$ID, rev(res2subtop15$ID))
p2 <- ggplot(res2subtop15, aes(x = ID, y = link))+
  geom_bar(stat = "identity", fill = "lightblue", color = "black", width = 0.7)+
  theme_classic2()+
  #scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()
p2
ggsave("../8.results/inuloc.multiloc.genes.pathway.pdf", p2, width = 5, height = 3)


########sample distribution of locations########
selorgas <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Saccharomyces cerevisiae", "Drosophila melanogaster")
usedlocs <- c("Nucleus", "Cytosol", "Membrane", "Mitochondrion", "Secreted")
colorinfo <- brewer.pal(5, "Set1")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- subloc[subloc$V2 %in% selorgas, ]

subloc$V2 <- factor(subloc$V2, selorgas)
p <- ggplot(subloc, aes(V2, fill=V3))+
  geom_bar(position = "stack")+
  scale_fill_manual(values = colorinfo)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p
pdf("../8.results/inuloc.sublocs.distribut.multiloc.organism.pdf", width = 5, height =4)
print(p)
dev.off()

########amino acid composition for all locs########
selorgas <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Saccharomyces cerevisiae", "Drosophila melanogaster")
usedlocs <- c("Nucleus", "Cytosol", "Membrane", "Mitochondrion", "Secreted")
AA_list_k=c("H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E")
colorinfo <- brewer.pal(5, "Set1")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")

suball <- subloc
subnuc <- subloc[subloc$V3 == "Nucleus", ]
subcyt <- subloc[subloc$V3 == "Cytosol", ]
submem <- subloc[subloc$V3 == "Membrane", ]
submit <- subloc[subloc$V3 == "Mitochondrion", ]
subsec <- subloc[subloc$V3 == "Secreted", ]

suballres <- as.data.frame(table(strsplit(paste(suball$V4, collapse = ""), "")[[1]]))
suballres$Rate <- suballres$Freq / sum(suballres$Freq)
subnucres <- as.data.frame(table(strsplit(paste(subnuc$V4, collapse = ""), "")[[1]]))
subnucres$Rate <- subnucres$Freq / sum(subnucres$Freq)
subcytres <- as.data.frame(table(strsplit(paste(subcyt$V4, collapse = ""), "")[[1]]))
subcytres$Rate <- subcytres$Freq / sum(subcytres$Freq)
submemres <- as.data.frame(table(strsplit(paste(submem$V4, collapse = ""), "")[[1]]))
submemres$Rate <- submemres$Freq / sum(submemres$Freq)
submitres <- as.data.frame(table(strsplit(paste(submit$V4, collapse = ""), "")[[1]]))
submitres$Rate <- submitres$Freq / sum(submitres$Freq)
subsecres <- as.data.frame(table(strsplit(paste(subsec$V4, collapse = ""), "")[[1]]))
subsecres$Rate <- subsecres$Freq / sum(subsecres$Freq)

res <- rbind(suballres, subnucres, subcytres, submemres, submitres, subsecres)
res$Class <- c(rep("BG", nrow(suballres)), rep("Nuc", nrow(subnucres)), rep("Cyt", nrow(subcytres)),
               rep("Mem", nrow(submemres)), rep("Mit", nrow(submitres)), rep("Sec", nrow(subsecres)))
res <- res[res$Var1 %in% AA_list_k, ]
res$Var1 <- factor(res$Var1, AA_list_k)
p6 <- ggplot(res, aes(x = Var1, y = Rate, group = Class))+
  geom_line(color = "gray50", alpha = 0.5)+
  geom_point(aes(shape=Class,colour=Var1), size = 5)+
  theme_classic2()
p6

pdf("../8.results/inuloc.sublocs.distribut.aa.freq.pdf", width = 10, height =7)
print(p6)
dev.off()

########TF mutation in NLS, substrate level, survival########
allnls <- fread("../0.dataset/nlsdb.nls.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsnls <- merge(hsdataset, allnls, by.x = "V1", by.y = "V1")
alltfs <- fread("../0.dataset/TF-Target-information.txt", header = T, stringsAsFactors = F, data.table = F)
hsnls <- hsnls[hsnls$V2.x %in% alltfs$TF, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)

Genes = c()
mutsamnum = c()
mutVsctrFC = c()
mutVsctrPva = c()
mutVsctrSurv = c()
mutVsctrSurvPva = c()
tumortype = c()
for (tfs in hsnls$V2.x) {
  nlsregion = hsnls[hsnls$V2.x == tfs, ]$V2.y
  regnum <- str2lst2(nlsregion, nchar(hsnls[hsnls$V2.x == tfs, ]$V3))
  subs <- alltfs[alltfs$TF == tfs, ]$target
  if (length(subs) < 3) {
    next
  }
  
  for (cancer in tumorlist$V1) {
    singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                     header = T, stringsAsFactors = F, data.table = F)
    expmtr <- fread(paste("/home/yukai/data/TCGA/htseq_fpkm/gene/TCGA-", cancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
    survfile <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", cancer, ".survival.tsv", sep = ""), 
                      header = T, stringsAsFactors = F, data.table = F)
    survfile$OS.days <- survfile$OS.time
    survfile$OS.status <- survfile$OS
    rownames(expmtr) <- expmtr$Ensembl_ID
    expmtr <- expmtr[, -1]
    expmtr <- expmtr[as.numeric(substr(names(expmtr), 14, 15)) < 10, ]
    
    singmut <- singmut[singmut$gene == tfs, ]
    mutsam <- singmut[as.numeric(substr(singmut$Amino_Acid_Change, 4, nchar(singmut$Amino_Acid_Change) -1)) %in% regnum, ]$Sample_ID
    
    if (length(mutsam) < 3) {
      Genes = c(Genes, tfs)
      mutsamnum = c(mutsamnum, length(mutsam))
      mutVsctrFC = c(mutVsctrFC, NA)
      mutVsctrPva = c(mutVsctrPva, NA)
      mutVsctrSurv = c(mutVsctrSurv, NA)
      mutVsctrSurvPva = c(mutVsctrSurvPva, NA)
      tumortype = c(tumortype, cancer)
    }else{
      Genes = c(Genes, tfs)
      mutsamnum = c(mutsamnum, length(mutsam))
      tumortype = c(tumortype, cancer)
      selexp1 = expmtr[rownames(expmtr) %in% subs, colnames(expmtr) %in% mutsam]
      selexp2 = expmtr[rownames(expmtr) %in% subs, !(colnames(expmtr) %in% mutsam)]
      mutVsctrFC = c(mutVsctrFC, log2(median(rowMedians(as.matrix(selexp1), na.rm = T))/median(rowMedians(as.matrix(selexp2), na.rm = T))))
      mutVsctrPva = c(mutVsctrPva, t.test(rowMedians(as.matrix(selexp1), na.rm = T), rowMedians(as.matrix(selexp2), na.rm = T))$p.value)
      survfile$Class = ifelse(survfile$sample %in% mutsam, "0.mut", "1.non")
      if (nrow(survfile[survfile$Class == "0.mut", ]) < 2) {
        mutVsctrSurv = c(mutVsctrSurv, 1)
        mutVsctrSurvPva <- c(mutVsctrSurvPva, 1)
      }else{
        tmp <- summary(coxph((Surv(OS.days, OS.status)) ~ Class, data = survfile))
        mutVsctrSurv = c(mutVsctrSurv, tmp$coefficients[[1]])
        mutVsctrSurvPva <- c(mutVsctrSurvPva, tmp$logtest[[3]])
      }
    }
  }
  print(paste(tfs, " in ", cancer, " has finished !!!!!!!"))
}

res <- data.frame(Genes = Genes,
                  Num = mutsamnum,
                  FC = mutVsctrFC,
                  Pva = mutVsctrPva,
                  Surv = mutVsctrSurv,
                  SurvPva = mutVsctrSurvPva,
                  Cancer = tumortype,
                  stringsAsFactors = F)

write.table(res, "../8.results/inuloc.point.mut.nls.txt", row.names = F, col.names = T, quote = F, sep = "\t")
########TF mutation in NLS, substrate level, survival for plot########
res <- fread("../8.results/inuloc.point.mut.nls.txt", header = T, stringsAsFactors = F, data.table = F)
tsg <- fread("/home/yukai/projects/CBAmodel/ProNuclear/0.database/Human_TSGs.txt", header = T, stringsAsFactors = F, data.table = F)
res <- res[res$Genes %in% tsg$GeneSymbol, ]

gene = "FOXA2"
cancer = "UCEC"
##mut lollipop
laml1 = read.maf(maf = paste('/home/yukai/data/TCGA/mutation_maf/', cancer, '.Mutation_filter.txt', sep = ""))
lollipopPlot(laml1, gene = gene, AACol = 'Protein_Change', labelPos = 'all', repel = T)
pdf('../8.results/inuloc.point.mut.nls.selectgene.lollip.pdf', width = 8, height = 5)
lollipopPlot(laml1, gene = gene, AACol = 'Protein_Change', labelPos = 'all', repel = T)
dev.off()
##mut and sub expression]
oncogene <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
nlsregion = hsnls[hsnls$V2.x == gene, ]$V2.y
regnum <- str2lst2(nlsregion, nchar(hsnls[hsnls$V2.x == tfs, ]$V3))
subs <- alltfs[alltfs$TF == gene, ]$target
subs <- intersect(subs,oncogene$Gene)
singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                 header = T, stringsAsFactors = F, data.table = F)
expmtr <- fread(paste("/home/yukai/data/TCGA/htseq_fpkm/gene/TCGA-", cancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]
expmtr <- expmtr[as.numeric(substr(names(expmtr), 14, 15)) < 10, ]
singmut <- singmut[singmut$gene == gene, ]
mutsam <- singmut[as.numeric(substr(singmut$Amino_Acid_Change, 4, nchar(singmut$Amino_Acid_Change) -1)) %in% regnum, ]$Sample_ID
selexp1 = expmtr[rownames(expmtr) %in% subs, colnames(expmtr) %in% mutsam]
selexp2 = expmtr[rownames(expmtr) %in% subs, !(colnames(expmtr) %in% mutsam)]

subs <- c()
logfc <- c()
pva <- c()
for (i in 1:nrow(selexp1)) {
  subs <- c(subs, rownames(selexp1)[i])
  logfc <- c(logfc, log2(median(as.numeric(selexp1[i, ])) / median(as.numeric(selexp2[i, ]))))
  pva <- c(pva, t.test(as.numeric(selexp1[i, ]), as.numeric(selexp2[i, ]))$p.value)
}
resde <- data.frame(ID = subs,
                    logFC = logfc,
                    Pva = pva)

targets <- c("CCL7", "CDC6")

selexp1res <- selexp1[targets, ]
selexp1res$Gene <- rownames(selexp1res)
selexp1res <- melt(selexp1res, id.vars = "Gene")
selexp1res$Type <- "Mut"
selexp2res <- selexp2[targets, ]
selexp2res$Gene <- rownames(selexp2res)
selexp2res <- melt(selexp2res, id.vars = "Gene")
selexp2res$Type <- "Other"

expres <- rbind(selexp1res, selexp2res)
p1 <- ggplot(expres, aes(x = Type, y = value, fill = Type))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", method = "wilcox")+
  ylab("Substrate expression level")+
  theme_classic2()
q <- facet(p1, facet.by = "Gene", scales = "free_y")
q
pdf("../8.results/inuloc.point.mut.nls.selectgene.subexp.pdf", width = 7, height =4)
print(q)
dev.off()
##survival
nlsregion = hsnls[hsnls$V2.x == gene, ]$V2.y
regnum <- str2lst2(nlsregion, nchar(hsnls[hsnls$V2.x == tfs, ]$V3))
subs <- alltfs[alltfs$TF == gene, ]$target
singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                 header = T, stringsAsFactors = F, data.table = F)
expmtr <- fread(paste("/home/yukai/data/TCGA/htseq_fpkm/gene/TCGA-", cancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                header = T, stringsAsFactors = F, data.table = F)
survfile <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", cancer, ".survival.tsv", sep = ""), 
                  header = T, stringsAsFactors = F, data.table = F)
survfile$OS.days <- survfile$OS.time
survfile$OS.status <- survfile$OS
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]
expmtr <- expmtr[as.numeric(substr(names(expmtr), 14, 15)) < 10, ]
singmut <- singmut[singmut$gene == gene, ]
mutsam <- singmut[as.numeric(substr(singmut$Amino_Acid_Change, 4, nchar(singmut$Amino_Acid_Change) -1)) %in% regnum, ]$Sample_ID
survfile$Class <- ifelse(survfile$sample %in% mutsam, "0.Mut", "1.Other")
tmp <- summary(coxph((Surv(OS.days, OS.status)) ~ Class, data = survfile))
fit <- survfit(Surv(OS.days, OS.status) ~ Class, data = survfile)
pval <- tmp$logtest[[3]]
p <- ggsurvplot(fit,
                linetype = 1,
                #censor.shape=45,
                size = 1, # change line size
                #palette = c("#6bb82c", "#e62019"),# custom color palettes
                conf.int = TRUE, # Add confidence interval
                pval = paste('p = ', pval), # Add p-value
                #risk.table.col = "strata",# Risk table color by groups
                legend = "top",
                xlab = "Time (days)",
                ylab = "Probability of OS")
p
pdf("../8.results/inuloc.point.mut.nls.selectgene.surv.pdf", width = 5, height =4)
print(p)
dev.off()

########human dataset annotation########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc[subloc$V2 == "Homo sapiens", ]
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]

sigregions1 <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions2 <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions <- rbind(sigregions1, sigregions2)
sigregions3 <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions4 <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)

res <- data.frame(ID = c("Nucleus", "Mito", "Secreted"),
                  Ratio = c(length(unique(intersect(sigregions$V1, subloc[subloc$V3 == "Nucleus", ]$V1))) / length(subloc[subloc$V3 == "Nucleus", ]$V1),
                            length(unique(intersect(sigregions4$V1, subloc[subloc$V3 == "Mitochondrion", ]$V1))) / length(subloc[subloc$V3 == "Mitochondrion", ]$V1),
                            length(unique(intersect(sigregions3$V1, subloc[subloc$V3 == "Secreted", ]$V1))) / length(subloc[subloc$V3 == "Secreted", ]$V1)))

p4 <- ggplot(res, aes(x = ID, y = Ratio, fill = ID))+
  geom_bar(stat="identity", position=position_dodge())+
  #stat_compare_means(label = "p.signif")+
  geom_hline(yintercept = 0.2, lty = "dashed")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.25))+
  ylab("Validated Regions")+
  theme_classic2()
p4
pdf("../8.results/inuloc.sublocs.distribut.cover.pdf", width = 4, height =2.5)
print(p4)
dev.off()

#################################Model Features#################################
########data summary########
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 == "Nucleus", ]
99920 - 36016

69946 - 25212
19984 - 7203
9992 - 3601

########10-fold cvs for Cytosol########
##different parameters prepare table
color10 <- brewer.pal(10, "Spectral")
cv01 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv10.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv12.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv13.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv14.txt", header =T, stringsAsFactors = F, data.table = F)
cv06 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv15.txt", header =T, stringsAsFactors = F, data.table = F)
cv07 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv16.txt", header =T, stringsAsFactors = F, data.table = F)
cv08 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv17.txt", header =T, stringsAsFactors = F, data.table = F)
cv09 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv18.txt", header =T, stringsAsFactors = F, data.table = F)
cv10 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv19.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                  round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                  round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.cv.Cytosol.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
lines.roc(cv06$label, cv06$score, col = color10[6], percent=TRUE)
lines.roc(cv07$label, cv07$score, col = color10[7], percent=TRUE)
lines.roc(cv08$label, cv08$score, col = color10[8], percent=TRUE)
lines.roc(cv09$label, cv09$score, col = color10[9], percent=TRUE)
lines.roc(cv10$label, cv10$score, col = color10[10], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##train test validation
color10 <- c("red", "blue", "green")
trdata <- fread("../2.TrainValidTest/roc_data.Cytosol.train.txt", header =T, stringsAsFactors = F, data.table = F)
vadata <- fread("../2.TrainValidTest/roc_data.Cytosol.valid.txt", header =T, stringsAsFactors = F, data.table = F)
tedata <- fread("../2.TrainValidTest/roc_data.Cytosol.test.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("Train", "Valid", "Test"),
                c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                  round(auc(tedata$label, tedata$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.trvate.Cytosol.pdf", width = 4, height = 4)
plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
dev.off()


########10-fold cvs for Mitochondrion########
##different parameters prepare table
color10 <- brewer.pal(10, "Spectral")
cv01 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv10.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv12.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv13.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv14.txt", header =T, stringsAsFactors = F, data.table = F)
cv06 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv15.txt", header =T, stringsAsFactors = F, data.table = F)
cv07 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv16.txt", header =T, stringsAsFactors = F, data.table = F)
cv08 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv17.txt", header =T, stringsAsFactors = F, data.table = F)
cv09 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv18.txt", header =T, stringsAsFactors = F, data.table = F)
cv10 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv19.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                  round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                  round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.cv.Mitochondrion.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
lines.roc(cv06$label, cv06$score, col = color10[6], percent=TRUE)
lines.roc(cv07$label, cv07$score, col = color10[7], percent=TRUE)
lines.roc(cv08$label, cv08$score, col = color10[8], percent=TRUE)
lines.roc(cv09$label, cv09$score, col = color10[9], percent=TRUE)
lines.roc(cv10$label, cv10$score, col = color10[10], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##train test validation
color10 <- c("red", "blue", "green")
trdata <- fread("../2.TrainValidTest/roc_data.Mitochondrion.train.txt", header =T, stringsAsFactors = F, data.table = F)
vadata <- fread("../2.TrainValidTest/roc_data.Mitochondrion.valid.txt", header =T, stringsAsFactors = F, data.table = F)
tedata <- fread("../2.TrainValidTest/roc_data.Mitochondrion.test.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("Train", "Valid", "Test"),
                c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                  round(auc(tedata$label, tedata$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.trvate.Mitochondrion.pdf", width = 4, height = 4)
plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
dev.off()


########10-fold cvs for Membrane########
##different parameters prepare table
color10 <- brewer.pal(10, "Spectral")
cv01 <- fread("../2.TrainValidTest/roc_data.Membrane.cv10.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../2.TrainValidTest/roc_data.Membrane.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../2.TrainValidTest/roc_data.Membrane.cv12.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../2.TrainValidTest/roc_data.Membrane.cv13.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../2.TrainValidTest/roc_data.Membrane.cv14.txt", header =T, stringsAsFactors = F, data.table = F)
cv06 <- fread("../2.TrainValidTest/roc_data.Membrane.cv15.txt", header =T, stringsAsFactors = F, data.table = F)
cv07 <- fread("../2.TrainValidTest/roc_data.Membrane.cv16.txt", header =T, stringsAsFactors = F, data.table = F)
cv08 <- fread("../2.TrainValidTest/roc_data.Membrane.cv17.txt", header =T, stringsAsFactors = F, data.table = F)
cv09 <- fread("../2.TrainValidTest/roc_data.Membrane.cv18.txt", header =T, stringsAsFactors = F, data.table = F)
cv10 <- fread("../2.TrainValidTest/roc_data.Membrane.cv19.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                  round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                  round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.cv.Membrane.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
lines.roc(cv06$label, cv06$score, col = color10[6], percent=TRUE)
lines.roc(cv07$label, cv07$score, col = color10[7], percent=TRUE)
lines.roc(cv08$label, cv08$score, col = color10[8], percent=TRUE)
lines.roc(cv09$label, cv09$score, col = color10[9], percent=TRUE)
lines.roc(cv10$label, cv10$score, col = color10[10], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##train test validation
color10 <- c("red", "blue", "green")
trdata <- fread("../2.TrainValidTest/roc_data.Membrane.train.txt", header =T, stringsAsFactors = F, data.table = F)
vadata <- fread("../2.TrainValidTest/roc_data.Membrane.valid.txt", header =T, stringsAsFactors = F, data.table = F)
tedata <- fread("../2.TrainValidTest/roc_data.Membrane.test.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("Train", "Valid", "Test"),
                c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                  round(auc(tedata$label, tedata$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.trvate.Membrane.pdf", width = 4, height = 4)
plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
dev.off()


########10-fold cvs for Secreted########
##different parameters prepare table
color10 <- brewer.pal(10, "Spectral")
cv01 <- fread("../2.TrainValidTest/roc_data.Secreted.cv10.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../2.TrainValidTest/roc_data.Secreted.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../2.TrainValidTest/roc_data.Secreted.cv12.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../2.TrainValidTest/roc_data.Secreted.cv13.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../2.TrainValidTest/roc_data.Secreted.cv14.txt", header =T, stringsAsFactors = F, data.table = F)
cv06 <- fread("../2.TrainValidTest/roc_data.Secreted.cv15.txt", header =T, stringsAsFactors = F, data.table = F)
cv07 <- fread("../2.TrainValidTest/roc_data.Secreted.cv16.txt", header =T, stringsAsFactors = F, data.table = F)
cv08 <- fread("../2.TrainValidTest/roc_data.Secreted.cv17.txt", header =T, stringsAsFactors = F, data.table = F)
cv09 <- fread("../2.TrainValidTest/roc_data.Secreted.cv18.txt", header =T, stringsAsFactors = F, data.table = F)
cv10 <- fread("../2.TrainValidTest/roc_data.Secreted.cv19.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                  round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                  round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.cv.Secreted.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
lines.roc(cv06$label, cv06$score, col = color10[6], percent=TRUE)
lines.roc(cv07$label, cv07$score, col = color10[7], percent=TRUE)
lines.roc(cv08$label, cv08$score, col = color10[8], percent=TRUE)
lines.roc(cv09$label, cv09$score, col = color10[9], percent=TRUE)
lines.roc(cv10$label, cv10$score, col = color10[10], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()
##train test validation
color10 <- c("red", "blue", "green")
trdata <- fread("../2.TrainValidTest/roc_data.Secreted.train.txt", header =T, stringsAsFactors = F, data.table = F)
vadata <- fread("../2.TrainValidTest/roc_data.Secreted.valid.txt", header =T, stringsAsFactors = F, data.table = F)
tedata <- fread("../2.TrainValidTest/roc_data.Secreted.test.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("Train", "Valid", "Test"),
                c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                  round(auc(tedata$label, tedata$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.trvate.Secreted.pdf", width = 4, height = 4)
plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
dev.off()



########10-fold cvs for combined########
color10 <- brewer.pal(5, "Spectral")
cv01 <- fread("../2.TrainValidTest/roc_data.Nucleus.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv02 <- fread("../2.TrainValidTest/roc_data.Cytosol.cv16.txt", header =T, stringsAsFactors = F, data.table = F)
cv03 <- fread("../2.TrainValidTest/roc_data.Mitochondrion.cv13.txt", header =T, stringsAsFactors = F, data.table = F)
cv04 <- fread("../2.TrainValidTest/roc_data.Membrane.cv11.txt", header =T, stringsAsFactors = F, data.table = F)
cv05 <- fread("../2.TrainValidTest/roc_data.Secreted.cv10.txt", header =T, stringsAsFactors = F, data.table = F)
cvaucs <- paste(c("cv_Nuc", "cv_Cyt", "cv_Mit", "cv_Mem", "cv_Sec"),
                c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                  round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                  round(auc(cv05$label, cv05$score), 4)),
                sep = " = ")
pdf("../8.results/model.feature.performance.cv.Combined.pdf", width = 4, height = 4)
plot.roc(cv01$label, cv01$score, col = color10[1], percent=TRUE)
lines.roc(cv02$label, cv02$score, col = color10[2], percent=TRUE)
lines.roc(cv03$label, cv03$score, col = color10[3], percent=TRUE)
lines.roc(cv04$label, cv04$score, col = color10[4], percent=TRUE)
lines.roc(cv05$label, cv05$score, col = color10[5], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=0.5)
dev.off()

########select data for comparison########
n = 500
#Cytosol         MULocDeep
data <- fread("../0.data/1.dataset.subloc.Cytosol.txt.test", header = T, stringsAsFactors = F, data.table = F)
data_pos <- data[data$Class == 1, ]
data_neg <- data[data$Class == 0, ]

set.seed(123456)
sub_nuc1 = sample(1:nrow(data_pos), n, replace=F)
valid_nuc = data_pos[sub_nuc1, ]

sub_no1 = sample(1:nrow(data_neg), n, replace=F)
valid_no = data_neg[sub_no1, ]

write.table(rbind(valid_nuc, valid_no), "../0.data/1.dataset.subloc.Cytosol.txt.compare", row.names = F, 
            col.names = T, sep = "\t", quote = F)


#Mitochondrion            DeepMito, MULocDeep
data <- fread("../0.data/1.dataset.subloc.Mitochondrion.txt.test", header = T, stringsAsFactors = F, data.table = F)
data_pos <- data[data$Class == 1, ]
data_neg <- data[data$Class == 0, ]

set.seed(123456)
sub_nuc1 = sample(1:nrow(data_pos), n, replace=F)
valid_nuc = data_pos[sub_nuc1, ]

sub_no1 = sample(1:nrow(data_neg), n, replace=F)
valid_no = data_neg[sub_no1, ]

write.table(rbind(valid_nuc, valid_no), "../0.data/1.dataset.subloc.Mitochondrion.txt.compare", row.names = F, 
            col.names = T, sep = "\t", quote = F)


#Membrane             MULocDeep
data <- fread("../0.data/1.dataset.subloc.Membrane.txt.test", header = T, stringsAsFactors = F, data.table = F)
data_pos <- data[data$Class == 1, ]
data_neg <- data[data$Class == 0, ]

set.seed(123456)
sub_nuc1 = sample(1:nrow(data_pos), n, replace=F)
valid_nuc = data_pos[sub_nuc1, ]

sub_no1 = sample(1:nrow(data_neg), n, replace=F)
valid_no = data_neg[sub_no1, ]

write.table(rbind(valid_nuc, valid_no), "../0.data/1.dataset.subloc.Membrane.txt.compare", row.names = F, 
            col.names = T, sep = "\t", quote = F)


#Secreted             SignalP, MULocDeep
data <- fread("../0.data/1.dataset.subloc.Secreted.txt.test", header = T, stringsAsFactors = F, data.table = F)
data_pos <- data[data$Class == 1, ]
data_neg <- data[data$Class == 0, ]

set.seed(123456)
sub_nuc1 = sample(1:nrow(data_pos), n, replace=F)
valid_nuc = data_pos[sub_nuc1, ]

sub_no1 = sample(1:nrow(data_neg), n, replace=F)
valid_no = data_neg[sub_no1, ]

write.table(rbind(valid_nuc, valid_no), "../0.data/1.dataset.subloc.Secreted.txt.compare", row.names = F, 
            col.names = T, sep = "\t", quote = F)

########compare with existing tools for four locs########
library("scales")
library("magrittr")
library("tidyr")
#Cytosol         MULocDeep
sysucc <- fread("../3.existingtools/sel.protein.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)
mulocdeep <- fread("../3.existingtools/mu-loc_Cyto.txt", header = F, stringsAsFactors = F, data.table = F)
mulocdeep$Score <- as.numeric(gsub("Cytoplasm:", "", mulocdeep$V3))
mulocdeep$Name <- gsub(":", "", gsub(">", "", mulocdeep$V1))
mulocdeep <- merge(mulocdeep, sysucc, by.x = "Name", by.y = "V1")
color10 <- brewer.pal(4, "Spectral")
cv01 <- sysucc
cv02 <- mulocdeep
cvaucs <- paste(c("SYSUCC", "MULocDeep"),
                c(round(auc(cv01$V6, cv01$V7), 3), 
                  round(auc(cv02$V6.y, cv02$Score), 3)),
                sep = " = ")
pdf("../3.existingtools/model.compare.existing.tools.Cytosol.pdf", width = 4, height = 4)
plot.roc(cv01$V6, cv01$V7, col = color10[1], percent=TRUE)
lines.roc(cv02$V6.y, cv02$Score, col = color10[2], percent=TRUE)
legend("bottomright", legend=cvaucs, col=c(color10[1:2]), lwd=2, cex=0.5)
dev.off()

#Mitochondrion            DeepMito, MULocDeep
sysucc <- fread("../3.existingtools/sel.protein.score.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
mulocdeep <- fread("../3.existingtools/mu-loc_Mito.txt", header = F, stringsAsFactors = F, data.table = F)
mulocdeep$Score <- as.numeric(gsub("Mitochondrion:", "", mulocdeep$V5))
mulocdeep$Name <- gsub(":", "", gsub(">", "", mulocdeep$V1))
mulocdeep <- merge(mulocdeep, sysucc, by.x = "Name", by.y = "V1")
deepmito <- fread("../3.existingtools/DeepMito.txt", header = T, stringsAsFactors = F, data.table = F)
deepmito <- merge(deepmito, sysucc, by.x = "Accession", by.y = "V1")
color10 <- brewer.pal(4, "Spectral")
cv01 <- sysucc
cv02 <- mulocdeep
cv03 <- deepmito
cvaucs <- paste(c("SYSUCC", "MULocDeep", "DeepMito"),
                c(round(auc(cv01$V6, cv01$V7), 3), 
                  round(auc(cv02$V6.y, cv02$Score), 3),
                  round(auc(cv03$V6, cv03$Score), 3)),
                sep = " = ")
pdf("../3.existingtools/model.compare.existing.tools.Mitochondrion.pdf", width = 4, height = 4)
plot.roc(cv01$V6, cv01$V7, col = color10[1], percent=TRUE)
lines.roc(cv02$V6.y, cv02$Score, col = color10[2], percent=TRUE)
lines.roc(cv03$V6, cv03$Score, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=c(color10[1:3]), lwd=2, cex=0.5)
dev.off()

#Membrane             MULocDeep
sysucc <- fread("../3.existingtools/sel.protein.score.Membrane", header = F, stringsAsFactors = F, data.table = F)
mulocdeep <- fread("../3.existingtools/mu-loc_Membrane.txt", header = F, stringsAsFactors = F, data.table = F)
mulocdeep$Score <- as.numeric(gsub("Membrane:", "", mulocdeep$V6))
mulocdeep$Name <- gsub(":", "", gsub(">", "", mulocdeep$V1))
mulocdeep <- merge(mulocdeep, sysucc, by.x = "Name", by.y = "V1")
color10 <- brewer.pal(4, "Spectral")
cv01 <- sysucc
cv02 <- mulocdeep
cvaucs <- paste(c("SYSUCC", "MULocDeep"),
                c(round(auc(cv01$V6, cv01$V7), 3), 
                  round(auc(cv02$V6.y, cv02$Score), 3)),
                sep = " = ")
pdf("../3.existingtools/model.compare.existing.tools.Membrane.pdf", width = 4, height = 4)
plot.roc(cv01$V6, cv01$V7, col = color10[1], percent=TRUE)
lines.roc(cv02$V6.y, cv02$Score, col = color10[2], percent=TRUE)
legend("bottomright", legend=cvaucs, col=c(color10[1:2]), lwd=2, cex=0.5)
dev.off()

#Secreted             SignalP, MULocDeep
sysucc <- fread("../3.existingtools/sel.protein.score.Secreted", header = F, stringsAsFactors = F, data.table = F)
mulocdeep <- fread("../3.existingtools/mu-loc_Secreted.txt", header = F, stringsAsFactors = F, data.table = F)
mulocdeep$Score <- as.numeric(gsub("Secreted:", "", mulocdeep$V4))
mulocdeep$Name <- gsub(":", "", gsub(">", "", mulocdeep$V1))
mulocdeep <- merge(mulocdeep, sysucc, by.x = "Name", by.y = "V1")
signalip <- fread("../3.existingtools/SignalP_Secreted.txt", header = T, stringsAsFactors = F, data.table = F)
signalip <- merge(signalip, sysucc, by.x = "# ID", by.y = "V1")
color10 <- brewer.pal(4, "Spectral")
cv01 <- sysucc
cv02 <- mulocdeep
cv03 <- signalip
cvaucs <- paste(c("SYSUCC", "MULocDeep", "SignalIP"),
                c(round(auc(cv01$V6, cv01$V7), 3), 
                  round(auc(cv02$V6.y, cv02$Score), 3),
                  round(auc(cv03$V6, cv03$`SP(Sec/SPI)`), 3)),
                sep = " = ")
pdf("../3.existingtools/model.compare.existing.tools.Secreted.pdf", width = 4, height = 4)
plot.roc(cv01$V6, cv01$V7, col = color10[1], percent=TRUE)
lines.roc(cv02$V6.y, cv02$Score, col = color10[2], percent=TRUE)
lines.roc(cv03$V6, cv03$`SP(Sec/SPI)`, col = color10[3], percent=TRUE)
legend("bottomright", legend=cvaucs, col=c(color10[1:3]), lwd=2, cex=0.5)
dev.off()


########SLPT for each loc########
hsdata <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
##Nucleus, Cytosol, Membrane, Mitochondrion, Secreted
cpt_nuc <- fread("../7.eachpro/all.CPT.delta.scores.Nucleus", header = F, stringsAsFactors = F, data.table = F)
cpt_cyt <- fread("../7.eachpro/all.CPT.delta.scores.Cytosol", header = F, stringsAsFactors = F, data.table = F)
cpt_mem <- fread("../7.eachpro/all.CPT.delta.scores.Membrane", header = F, stringsAsFactors = F, data.table = F)
cpt_mit <- fread("../7.eachpro/all.CPT.delta.scores.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
cpt_sec <- fread("../7.eachpro/all.CPT.delta.scores.Secreted", header = F, stringsAsFactors = F, data.table = F)

cpt_nuc <- cpt_nuc[cpt_nuc$V1 %in% hsdata$V1, ]
cpt_cyt <- cpt_cyt[cpt_cyt$V1 %in% hsdata$V1, ]
cpt_mem <- cpt_mem[cpt_mem$V1 %in% hsdata$V1, ]
cpt_mit <- cpt_mit[cpt_mit$V1 %in% hsdata$V1, ]
cpt_sec <- cpt_sec[cpt_sec$V1 %in% hsdata$V1, ]

write.table(cpt_nuc, "../7.eachpro/hs.CPT.delta.scores.Nucleus", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_cyt, "../7.eachpro/hs.CPT.delta.scores.Cytosol", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_mem, "../7.eachpro/hs.CPT.delta.scores.Membrane", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_mit, "../7.eachpro/hs.CPT.delta.scores.Mitochondrion", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_sec, "../7.eachpro/hs.CPT.delta.scores.Secreted", row.names = F, col.names = T, sep = "\t", quote = F)


cpt_nuc <- fread("../7.eachpro/all.CPT.delta.scores.regions.Nucleus", header = F, stringsAsFactors = F, data.table = F)
cpt_cyt <- fread("../7.eachpro/all.CPT.delta.scores.regions.Cytosol", header = F, stringsAsFactors = F, data.table = F)
cpt_mem <- fread("../7.eachpro/all.CPT.delta.scores.regions.Membrane", header = F, stringsAsFactors = F, data.table = F)
cpt_mit <- fread("../7.eachpro/all.CPT.delta.scores.regions.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
cpt_sec <- fread("../7.eachpro/all.CPT.delta.scores.regions.Secreted", header = F, stringsAsFactors = F, data.table = F)

cpt_nuc <- cpt_nuc[cpt_nuc$V1 %in% hsdata$V1, ]
cpt_cyt <- cpt_cyt[cpt_cyt$V1 %in% hsdata$V1, ]
cpt_mem <- cpt_mem[cpt_mem$V1 %in% hsdata$V1, ]
cpt_mit <- cpt_mit[cpt_mit$V1 %in% hsdata$V1, ]
cpt_sec <- cpt_sec[cpt_sec$V1 %in% hsdata$V1, ]

write.table(cpt_nuc, "../7.eachpro/hs.CPT.delta.scores.regions.Nucleus", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_cyt, "../7.eachpro/hs.CPT.delta.scores.regions.Cytosol", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_mem, "../7.eachpro/hs.CPT.delta.scores.regions.Membrane", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_mit, "../7.eachpro/hs.CPT.delta.scores.regions.Mitochondrion", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cpt_sec, "../7.eachpro/hs.CPT.delta.scores.regions.Secreted", row.names = F, col.names = T, sep = "\t", quote = F)

########SLPT trajectory for 600(200 for each) selected proteins########
hsdata <- fread("../4.panMut/hs.protein.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)
posdata <- fread("../0.data/0.2Cytosol_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
negdata <- fread("../0.data/0.2Cytosol_negative_data.txt", header = F, stringsAsFactors = F, data.table = F)
nondata <- setdiff(hsdata$V1, union(posdata$V1, negdata$V1))
res <- data.frame(ID = c(nondata, posdata$V1, negdata$V1),
                  LocInfo = c(rep("NoLoc", length(nondata)), rep("Nucle", length(posdata$V1)), rep("NonNucle", length(negdata$V1))))
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Cytosol", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
sam1 = nondata
sam2 = intersect(posdata$V1, unique(cptdata$UniID))
sam3 = intersect(negdata$V1, unique(cptdata$UniID))
set.seed(1234)
loc1 = sample(1:length(sam1), 500, replace=F)
loc2 = sample(1:length(sam2), 500, replace=F)
loc3 = sample(1:length(sam3), 500, replace=F)
uni1 = sam1[loc1]
uni2 = sam2[loc2]
uni3 = sam3[loc3]
cptdata4plot = cptdata[cptdata$UniID %in% c(uni1, uni2, uni3), ]
cptdata4plot <- merge(cptdata4plot, res, by.x = "UniID", by.y = "ID")
cptdata4plot$Strunc <- as.numeric(cptdata4plot$Strunc)
cptdata4plot$Pos <- as.numeric(cptdata4plot$Pos)
##allprotein
rolling_median <- function(formula, data, n_roll = 11, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}
predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}
p1 <- ggplot() +
  geom_line(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo, group = UniID), alpha = 0.1) +
  geom_smooth(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo), formula = y ~ x, method = "rolling_median", se = FALSE)+
  labs(title="Cytosol Location Probability Trajectory of all protein")+
  xlab("Protein sequence length")+
  ylab("Cytosol Location Probability")+
  scale_color_manual(values = c("green", "blue", "red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../8.results/NLPT.overview.Cytosol.pdf", p1, width = 10, height = 4)

hsdata <- fread("../4.panMut/hs.protein.score.Membrane", header = F, stringsAsFactors = F, data.table = F)
posdata <- fread("../0.data/0.2Membrane_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
negdata <- fread("../0.data/0.2Membrane_negative_data.txt", header = F, stringsAsFactors = F, data.table = F)
nondata <- setdiff(hsdata$V1, union(posdata$V1, negdata$V1))
res <- data.frame(ID = c(nondata, posdata$V1, negdata$V1),
                  LocInfo = c(rep("NoLoc", length(nondata)), rep("Nucle", length(posdata$V1)), rep("NonNucle", length(negdata$V1))))
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Membrane", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
sam1 = nondata
sam2 = intersect(posdata$V1, unique(cptdata$UniID))
sam3 = intersect(negdata$V1, unique(cptdata$UniID))
set.seed(12345)
loc1 = sample(1:length(sam1), 500, replace=F)
loc2 = sample(1:length(sam2), 500, replace=F)
loc3 = sample(1:length(sam3), 500, replace=F)
uni1 = sam1[loc1]
uni2 = sam2[loc2]
uni3 = sam3[loc3]
cptdata4plot = cptdata[cptdata$UniID %in% c(uni1, uni2, uni3), ]
cptdata4plot <- merge(cptdata4plot, res, by.x = "UniID", by.y = "ID")
cptdata4plot$Strunc <- as.numeric(cptdata4plot$Strunc)
cptdata4plot$Pos <- as.numeric(cptdata4plot$Pos)
##allprotein
rolling_median <- function(formula, data, n_roll = 11, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}
predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}
p1 <- ggplot() +
  geom_line(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo, group = UniID), alpha = 0.1) +
  geom_smooth(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo), formula = y ~ x, method = "rolling_median", se = FALSE)+
  labs(title="Membrane Location Probability Trajectory of all protein")+
  xlab("Protein sequence length")+
  ylab("Membrane Location Probability")+
  scale_color_manual(values = c("green", "blue", "red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../8.results/NLPT.overview.Membrane.pdf", p1, width = 10, height = 4)


hsdata <- fread("../4.panMut/hs.protein.score.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
posdata <- fread("../0.data/0.2Mitochondrion_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
negdata <- fread("../0.data/0.2Mitochondrion_negative_data.txt", header = F, stringsAsFactors = F, data.table = F)
nondata <- setdiff(hsdata$V1, union(posdata$V1, negdata$V1))
res <- data.frame(ID = c(nondata, posdata$V1, negdata$V1),
                  LocInfo = c(rep("NoLoc", length(nondata)), rep("Nucle", length(posdata$V1)), rep("NonNucle", length(negdata$V1))))
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Mitochondrion", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
sam1 = nondata
sam2 = intersect(posdata$V1, unique(cptdata$UniID))
sam3 = intersect(negdata$V1, unique(cptdata$UniID))
set.seed(12345678)
loc1 = sample(1:length(sam1), 500, replace=F)
loc2 = sample(1:length(sam2), 500, replace=F)
loc3 = sample(1:length(sam3), 500, replace=F)
uni1 = sam1[loc1]
uni2 = sam2[loc2]
uni3 = sam3[loc3]
cptdata4plot = cptdata[cptdata$UniID %in% c(uni1, uni2, uni3), ]
cptdata4plot <- merge(cptdata4plot, res, by.x = "UniID", by.y = "ID")
cptdata4plot$Strunc <- as.numeric(cptdata4plot$Strunc)
cptdata4plot$Pos <- as.numeric(cptdata4plot$Pos)
##allprotein
rolling_median <- function(formula, data, n_roll = 11, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}
predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}
p1 <- ggplot() +
  geom_line(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo, group = UniID), alpha = 0.1) +
  geom_smooth(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo), formula = y ~ x, method = "rolling_median", se = FALSE)+
  labs(title="Mitochondrion Location Probability Trajectory of all protein")+
  xlab("Protein sequence length")+
  ylab("Mitochondrion Location Probability")+
  scale_color_manual(values = c("green", "blue", "red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../8.results/NLPT.overview.Mitochondrion.pdf", p1, width = 10, height = 4)


hsdata <- fread("../4.panMut/hs.protein.score.Secreted", header = F, stringsAsFactors = F, data.table = F)
posdata <- fread("../0.data/0.2Secreted_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
negdata <- fread("../0.data/0.2Secreted_negative_data.txt", header = F, stringsAsFactors = F, data.table = F)
nondata <- setdiff(hsdata$V1, union(posdata$V1, negdata$V1))
res <- data.frame(ID = c(nondata, posdata$V1, negdata$V1),
                  LocInfo = c(rep("NoLoc", length(nondata)), rep("Nucle", length(posdata$V1)), rep("NonNucle", length(negdata$V1))))
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Secreted", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
sam1 = nondata
sam2 = intersect(posdata$V1, unique(cptdata$UniID))
sam3 = intersect(negdata$V1, unique(cptdata$UniID))
set.seed(1234)
loc1 = sample(1:length(sam1), 500, replace=F)
loc2 = sample(1:length(sam2), 500, replace=F)
loc3 = sample(1:length(sam3), 500, replace=F)
uni1 = sam1[loc1]
uni2 = sam2[loc2]
uni3 = sam3[loc3]
cptdata4plot = cptdata[cptdata$UniID %in% c(uni1, uni2, uni3), ]
cptdata4plot <- merge(cptdata4plot, res, by.x = "UniID", by.y = "ID")
cptdata4plot$Strunc <- as.numeric(cptdata4plot$Strunc)
cptdata4plot$Pos <- as.numeric(cptdata4plot$Pos)
##allprotein
rolling_median <- function(formula, data, n_roll = 11, ...) {
  x <- data$x[order(data$x)]
  y <- data$y[order(data$x)]
  y <- zoo::rollmedian(y, n_roll, na.pad = TRUE)
  structure(list(x = x, y = y, f = approxfun(x, y)), class = "rollmed")
}
predict.rollmed <- function(mod, newdata, ...) {
  setNames(mod$f(newdata$x), newdata$x)
}
p1 <- ggplot() +
  geom_line(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo, group = UniID), alpha = 0.1) +
  geom_smooth(data = cptdata4plot,aes(x = Pos,y = Strunc, color = LocInfo), formula = y ~ x, method = "rolling_median", se = FALSE)+
  labs(title="Secreted Location Probability Trajectory of all protein")+
  xlab("Protein sequence length")+
  ylab("Secreted Location Probability")+
  scale_color_manual(values = c("green", "blue", "red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
ggsave("../8.results/NLPT.overview.Secreted.pdf", p1, width = 10, height = 4)

########signal peptide delta scores distribution########
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Mitochondrion", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.dataset/mito.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
cptdata <- cptdata[cptdata$UniID %in% nlss$V1, ]
nlss$V3 <- as.numeric(gsub("\\..*", "", nlss$V2))
nlss$V4 <- as.numeric(gsub(".*\\.", "", nlss$V2))
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & as.numeric(cptdata$Pos) >= nlss[i, ]$V3 & as.numeric(cptdata$Pos) <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res
cptnuc <- cptdata
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = "NLS"
cptnuc4random$Type = "Random"
finalres <- rbind(resnuc, cptnuc4random)
finalres$Delta <- as.numeric(finalres$Delta)
p1 <- ggplot(finalres, aes(x=Type, y=Delta, color = Type))+
  geom_boxplot()+
  xlab("")+
  ylab("Delta Scores")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
ggsave("../8.results/NLPT.Mitochondrion.deltascore.pdf", p1, height = 4, width = 4)


cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Secreted", header =F, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.dataset/sp.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
cptdata <- cptdata[cptdata$UniID %in% nlss$V1, ]
nlss$V3 <- as.numeric(gsub("\\..*", "", nlss$V2))
nlss$V4 <- as.numeric(gsub(".*\\.", "", nlss$V2))
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & as.numeric(cptdata$Pos) >= nlss[i, ]$V3 & as.numeric(cptdata$Pos) <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res
cptnuc <- cptdata
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = "NLS"
cptnuc4random$Type = "Random"
finalres <- rbind(resnuc, cptnuc4random)
finalres$Delta <- as.numeric(finalres$Delta)
p1 <- ggplot(finalres, aes(x=Type, y=Delta, color = Type))+
  geom_boxplot()+
  xlab("")+
  ylab("Delta Scores")+
  stat_compare_means(method = "t.test")+
  theme_classic2()
p1
ggsave("../8.results/NLPT.Secreted.deltascore.pdf", p1, height = 4, width = 4)

########Sigdelta region retrive########
##Nucleus, Cytosol, Membrane, Mitochondrion, Secreted
cptdata1 <- fread("../7.eachpro/hs.CPT.delta.scores.Cytosol", header =T, stringsAsFactors = F, data.table = F)
cptdata1 <- cptdata1[as.numeric(cptdata1$V8) > 0.1, ]
write.table(cptdata1, "../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol", row.names = F, col.names = T, sep = "\t", quote = F)
cptdata1 <- fread("../7.eachpro/hs.CPT.delta.scores.Membrane", header =T, stringsAsFactors = F, data.table = F)
cptdata1 <- cptdata1[as.numeric(cptdata1$V8) > 0.1, ]
write.table(cptdata1, "../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane", row.names = F, col.names = T, sep = "\t", quote = F)
cptdata1 <- fread("../7.eachpro/hs.CPT.delta.scores.Mitochondrion", header =T, stringsAsFactors = F, data.table = F)
cptdata1 <- cptdata1[as.numeric(cptdata1$V8) > 0.1, ]
write.table(cptdata1, "../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion", row.names = F, col.names = T, sep = "\t", quote = F)
cptdata1 <- fread("../7.eachpro/hs.CPT.delta.scores.Secreted", header =T, stringsAsFactors = F, data.table = F)
cptdata1 <- cptdata1[as.numeric(cptdata1$V8) > 0.1, ]
write.table(cptdata1, "../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted", row.names = F, col.names = T, sep = "\t", quote = F)

########Truncated ternimal scores########
alllocs <- c("Cytosol", "Membrane", "Mitochondrion", "Secreted")
for (locs in alllocs) {
  hspos <- fread(paste("../0.data/0.2", locs, "_positive_data.txt", sep = ""), header = F, stringsAsFactors = F, data.table = F)
  hsdataset <- fread(paste("../4.panMut/hs.protein.score.", locs, sep = ""), header =F, stringsAsFactors = F, data.table = F)
  names(hsdataset) <- c("UniID", "Gene", "Seq", "Score")
  hstruncated <- fread(paste("../6.panStrunc/human.truncated.dataset.scores.", locs, sep = ""), header =T, stringsAsFactors = F)
  
  hsres <- hsdataset
  hsres$TruncatedTerminal <- "WT"
  hsres <- hsres[hsres$UniID %in% hspos$V1, ]
  res <- rbind(hsres[, c(1, 5, 4)], hstruncated[, c(1, 6, 7)])
  my_comparisons <- list(c("C", "WT"), c("M", "WT"), c("N", "WT"))
  p2 <- ggplot(res, aes(x=TruncatedTerminal, y=Score, color = TruncatedTerminal))+
    geom_boxplot()+
    xlab("")+
    ylab("Predicted Scores")+
    stat_compare_means(method = "t.test", comparisons = my_comparisons)+
    theme_classic2()
  p2
  ggsave(paste("../8.results/human.dataset.ternimal.truncated.", locs, ".score.pdf", sep = ""), p2, width = 6, height = 4)
  
}

########Truncated signal peptide and sig delta regions########
#Cytosol
hspos <- fread("../0.data/0.2Cytosol_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Cytosol", header =F, stringsAsFactors = F, data.table = F)
names(hsdataset) <- c("UniID", "Gene", "Seq", "Score")

nlsdata <- fread("../6.panStrunc/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres <- nlsdata[, c(1, 4, 8, 9)]
names(nlsres) <- c("UniID", "WT", "Truncted", "Random")
nlsres <- melt(nlsres, id.vars = "UniID")
p2 <- ggplot(nlsres, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=UniID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.05)+
  ggtitle("NLS region truncated")+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means()+
  theme_classic2()
p2
pdf("../8.results/human.dataset.Cytosol.truncated.pdf", width = 5, height = 4)
print(p2)
dev.off()

#Membrane
hspos <- fread("../0.data/0.2Membrane_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Membrane", header =F, stringsAsFactors = F, data.table = F)
names(hsdataset) <- c("UniID", "Gene", "Seq", "Score")

nlsdata <- fread("../6.panStrunc/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres <- nlsdata[, c(1, 4, 8, 9)]
names(nlsres) <- c("UniID", "WT", "Truncted", "Random")
nlsres <- melt(nlsres, id.vars = "UniID")
p2 <- ggplot(nlsres, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=UniID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.05)+
  ggtitle("NLS region truncated")+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means()+
  theme_classic2()
p2
pdf("../8.results/human.dataset.Membrane.truncated.pdf", width = 5, height = 4)
print(p2)
dev.off()

#Mitochondrion
nlsdata1 <- fread("../6.panStrunc/mito.human.pos.txt.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres1 <- nlsdata1[, c(1, 4, 8, 9)]
names(nlsres1) <- c("UniID.nls", "WT.nls", "Truncted.nls", "Random.nls")
nlsdata2 <- fread("../6.panStrunc/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres2 <- nlsdata2[, c(1, 4, 8, 9)]
names(nlsres2) <- c("UniID.yu", "WT.yu", "Truncted.yu", "Random.yu")
nlsresmerge <- merge(nlsres1, nlsres2, by.x = "UniID.nls", by.y = "UniID.yu")
nlsresmerge <- nlsresmerge[, c(1,2,3,6,4)]
names(nlsresmerge) <- c("ID", "WT", "NLS.Trunc", "Yu.Trunc", "Random")
nlsmelt <- melt(nlsresmerge, id.vars = "ID")
my_comparisons <- list( c("WT", "NLS.Trunc"), c("NLS.Trunc", "Yu.Trunc"),
                        c("WT", "Yu.Trunc"), c("WT", "Random"))
p3 <- ggplot(nlsmelt, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=ID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.1)+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic2()
p3
pdf("../8.results/human.dataset.Mitochondrion.truncated.pdf", width = 6, height = 4)
print(p3)
dev.off()

#Secreted
nlsdata1 <- fread("../6.panStrunc/sp.human.pos.txt.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres1 <- nlsdata1[, c(1, 4, 8, 9)]
names(nlsres1) <- c("UniID.nls", "WT.nls", "Truncted.nls", "Random.nls")
nlsdata2 <- fread("../6.panStrunc/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.4predict.results", header =T, stringsAsFactors = F, data.table = F)
nlsres2 <- nlsdata2[, c(1, 4, 8, 9)]
names(nlsres2) <- c("UniID.yu", "WT.yu", "Truncted.yu", "Random.yu")
nlsresmerge <- merge(nlsres1, nlsres2, by.x = "UniID.nls", by.y = "UniID.yu")
nlsresmerge <- nlsresmerge[, c(1,2,3,6,4)]
names(nlsresmerge) <- c("ID", "WT", "NLS.Trunc", "Yu.Trunc", "Random")
nlsmelt <- melt(nlsresmerge, id.vars = "ID")
my_comparisons <- list( c("WT", "NLS.Trunc"), c("NLS.Trunc", "Yu.Trunc"),
                        c("WT", "Yu.Trunc"), c("WT", "Random"))
p3 <- ggplot(nlsmelt, aes(x=variable, y=value, color = variable))+
  geom_boxplot()+
  geom_line(aes(group=ID), linetype="dashed", col="gray50", lwd=0.05, alpha = 0.1)+
  xlab("")+
  ylab("Predicted Scores")+
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic2()
p3
pdf("../8.results/human.dataset.Secreted.truncated.pdf", width = 6, height = 4)
print(p3)
dev.off()

########high delta terminal########
alllocs <- c("Cytosol", "Membrane", "Mitochondrion", "Secreted")
for (locs in alllocs) {
  hspos <- fread(paste("../0.data/0.2", locs, "_positive_data.txt", sep = ""), header = F, stringsAsFactors = F, data.table = F)
  cptdata <- fread(paste("../7.eachpro/hs.CPT.delta.scores.", locs, sep = ""), header =T, stringsAsFactors = F, data.table = F)
  names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
  cptdata <- cptdata[cptdata$Delta > 0.1, ]
  cptdata <- cptdata[cptdata$UniID %in% hspos$V1, ]
  uniid <- c()
  pos <- c()
  score <- c()
  highest <- c()
  for (ids in unique(cptdata$UniID)) {
    tmpdata <- cptdata[cptdata$UniID == ids, ]
    highest <- c(highest, max(tmpdata$Pos))
    tmpdata <- tmpdata[tmpdata$Pos == min(tmpdata$Pos), ]
    uniid <- c(uniid, tmpdata$UniID)
    pos <- c(pos, tmpdata$Pos)
    score <- c(score, tmpdata$Delta)
  }
  res <- data.frame(ID = uniid,
                    Pos = pos,
                    Score = score,
                    Length = highest)
  res$rela <- res$Pos/res$Length
  p <- ggplot(res, aes(x=rela))+
    geom_histogram(bins = 30,alpha=0.5,colour="black",size=0.25)+#, aes(fill = ..count..) )
    xlab("Relative Position to N terminal")+
    theme_bw()
  p
  ggsave(paste("../8.results/human.dataset.", locs, ".relaNterm.pdf", sep = ""), p, width = 5, height = 4)
  
}

########AA composition########
nlsseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.nls.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
nesseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.nes.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
cytseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
memseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
mitseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
secseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
bgsseq <- read.delim2("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.background.fasta",
                      col.names = "NlsSeqs", comment.char = ">")
seq2count <- function(seqs, group){
  akl <- nchar(as.character(seqs$NlsSeqs))
  akl_sum <- sum(akl, na.rm = T)
  AA_list_k=c("H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E")
  knum_sum <- c()
  for (AA in AA_list_k) {
    knum <- str_count(seqs$NlsSeqs, AA)
    knum_sum <- c(knum_sum, sum(knum, na.rm = T))
  }
  res <- data.frame(ID = AA_list_k,
                    Sum = knum_sum,
                    AllSum = akl_sum,
                    Ratio = knum_sum/akl_sum,
                    Class = group)
}

res <- rbind(seq2count(nlsseq, "NLS"), seq2count(nesseq, "NES"), seq2count(cytseq, "Cytosol"), 
             seq2count(memseq, "Membrane"), seq2count(mitseq, "Mitochondrion"), seq2count(secseq, "Secreted"),
             seq2count(bgsseq, "BGs"))
res$ID <- factor(res$ID, c("H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E"))
res$Class <- factor(res$Class, c("NLS", "NES", "Cytosol", "Membrane", "Mitochondrion", "Secreted", "BGs"))
p <- ggplot(res, aes(x = Class, y = ID, fill = Ratio))+
  geom_tile(size = 0.5, na.rm = T, color = "gray50")+
  scale_fill_gradient2(low="#B6C9F0", high="#E93B81",mid = "#FFE5E2", midpoint = 0.1)+
  coord_flip()+
  theme_classic2()
p
ggsave("../8.results/human.dataset.regions.AAcount.pdf", p, width = 7, height = 3)


########annotation rates########
usedlocs <- c("Nucleus", "Cytosol", "Mitochondrion", "Membrane", "Secreted")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc[subloc$V2 == "Homo sapiens", ]
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
sigregions1 <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions2 <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions <- rbind(sigregions1, sigregions2)
sigregions3 <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions4 <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
res1 <- data.frame(ID = c("Nucleus", "Mito", "Secreted", "Cytosol", "Membrane"),
                  Class = "ExpValid",
                  Ratio = c(length(unique(intersect(sigregions$V1, subloc[subloc$V3 == "Nucleus", ]$V1))) / length(subloc[subloc$V3 == "Nucleus", ]$V1),
                            length(unique(intersect(sigregions4$V1, subloc[subloc$V3 == "Mitochondrion", ]$V1))) / length(subloc[subloc$V3 == "Mitochondrion", ]$V1),
                            length(unique(intersect(sigregions3$V1, subloc[subloc$V3 == "Secreted", ]$V1))) / length(subloc[subloc$V3 == "Secreted", ]$V1),
                            0, 0))


preregions1 <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Nucleus.regions", header = F, stringsAsFactors = F, data.table = F)
preregions2 <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions", header = F, stringsAsFactors = F, data.table = F)
preregions3 <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions", header = F, stringsAsFactors = F, data.table = F)
preregions4 <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions", header = F, stringsAsFactors = F, data.table = F)
preregions5 <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions", header = F, stringsAsFactors = F, data.table = F)

res2 <- data.frame(ID = c("Nucleus", "Cytosol", "Mito", "Membrane", "Secreted"),
                   Class = "Predict",
                   Ratio = c(length(unique(intersect(preregions1$V1, subloc[subloc$V3 == "Nucleus", ]$V1))) / length(subloc[subloc$V3 == "Nucleus", ]$V1),
                             length(unique(intersect(preregions2$V1, subloc[subloc$V3 == "Cytosol", ]$V1))) / length(subloc[subloc$V3 == "Cytosol", ]$V1),
                             length(unique(intersect(preregions3$V1, subloc[subloc$V3 == "Mitochondrion", ]$V1))) / length(subloc[subloc$V3 == "Mitochondrion", ]$V1),
                             length(unique(intersect(preregions4$V1, subloc[subloc$V3 == "Membrane", ]$V1))) / length(subloc[subloc$V3 == "Membrane", ]$V1),
                             length(unique(intersect(preregions5$V1, subloc[subloc$V3 == "Secreted", ]$V1))) / length(subloc[subloc$V3 == "Secreted", ]$V1)))

res <- rbind(res1, res2)
res$ID <- factor(res$ID, c("Mito", "Nucleus", "Secreted", "Cytosol", "Membrane"))
p4 <- ggplot(res, aes(x = ID, y = Ratio, fill = Class))+
  geom_bar(stat="identity", position=position_dodge())+
  #stat_compare_means(label = "p.signif")+
  geom_hline(yintercept = median(res2$Ratio), lty = "dashed")+
  scale_y_continuous(expand = c(0, 0))+
  ylab("Annotated Regions")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p4
pdf("../8.results/human.dataset.annotation.pdf", width = 4, height =2.5)
print(p4)
dev.off()

########known regions ROCs########
cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Nucleus", header =T, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.dataset/mito.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
cptdata <- cptdata[cptdata$UniID %in% nlss$V1, ]
nlss$V3 <- as.numeric(gsub("\\..*", "", nlss$V2))
nlss$V4 <- as.numeric(gsub(".*\\.", "", nlss$V2))
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & as.numeric(cptdata$Pos) >= nlss[i, ]$V3 & as.numeric(cptdata$Pos) <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res
cptnuc <- cptdata
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = 1
cptnuc4random$Type = 0
finalres1 <- rbind(resnuc, cptnuc4random)


cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Mitochondrion", header =T, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.dataset/mito.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
cptdata <- cptdata[cptdata$UniID %in% nlss$V1, ]
nlss$V3 <- as.numeric(gsub("\\..*", "", nlss$V2))
nlss$V4 <- as.numeric(gsub(".*\\.", "", nlss$V2))
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & as.numeric(cptdata$Pos) >= nlss[i, ]$V3 & as.numeric(cptdata$Pos) <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res
cptnuc <- cptdata
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = 1
cptnuc4random$Type = 0
finalres2 <- rbind(resnuc, cptnuc4random)

cptdata <- fread("../7.eachpro/hs.CPT.delta.scores.Secreted", header =T, stringsAsFactors = F, data.table = F)
names(cptdata) <- c("UniID", "UniName", "GeneName", "Species", "Pos", "Score", "Strunc", "Delta")
nlss <- fread("../0.dataset/sp.human.pos.txt", header =F, stringsAsFactors = F, data.table = F)
cptdata <- cptdata[cptdata$UniID %in% nlss$V1, ]
nlss$V3 <- as.numeric(gsub("\\..*", "", nlss$V2))
nlss$V4 <- as.numeric(gsub(".*\\.", "", nlss$V2))
res <- cptdata[1,]
res <- res[-1, ]
for (i in c(1:nrow(nlss))) {
  tmpres <- cptdata[cptdata$UniID == nlss[i, ]$V1 & as.numeric(cptdata$Pos) >= nlss[i, ]$V3 & as.numeric(cptdata$Pos) <= nlss[i, ]$V4, ]
  res <- rbind(res, tmpres)
}
resnuc <- res
cptnuc <- cptdata
set.seed(123456)
nums = sample(1:nrow(cptnuc), nrow(resnuc), replace=F)
cptnuc4random <- cptnuc[nums, ]
resnuc$Type = 1
cptnuc4random$Type = 0
finalres3 <- rbind(resnuc, cptnuc4random)

color10 <- brewer.pal(10, "Spectral")
cvaucs <- paste(c("Nucleus", "Mito", "Secreted"),
                c(round(auc(finalres1$Type, finalres1$Delta), 4), round(auc(finalres2$Type, finalres2$Delta), 4),
                  round(auc(finalres3$Type, finalres3$Delta), 4)),
                sep = " = ")
pdf("../8.results/NLPT.deltascore.valid.regions.roc.pdf", width = 4, height = 4)
plot.roc(finalres1$Type, finalres1$Delta, col = color10[1], percent=TRUE)
lines.roc(finalres2$Type, finalres2$Delta, col = color10[3], percent=TRUE)
lines.roc(finalres3$Type, finalres3$Delta, col = color10[5], percent=TRUE)
legend("bottomright", legend=cvaucs, col=color10[c(1,3,5)], lwd=1, cex=0.5)
dev.off()

########motif enrichments########
#nothing, commond
########pSTY functional score########
res <- fread("../7.eachpro/pSTY.NLS.position.csv", header = T, stringsAsFactors = F, data.table = F, sep = ",")
res <- res[res$Fall_on != "NaN", ]
p1 <- ggplot(res, aes(x = Fall_on, y = Functional_score, fill = Fall_on))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("NLS")+
  ylab("Functional Score")+
  theme_classic2()
p1

res <- fread("../7.eachpro/pSTY.Cytosol.position.csv", header = T, stringsAsFactors = F, data.table = F, sep = ",")
res <- res[res$Fall_on != "NaN", ]
p2 <- ggplot(res, aes(x = Fall_on, y = Functional_score, fill = Fall_on))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("Cytosol")+
  ylab("Functional Score")+
  theme_classic2()
p2

res <- fread("../7.eachpro/pSTY.Membrane.position.csv", header = T, stringsAsFactors = F, data.table = F, sep = ",")
res <- res[res$Fall_on != "NaN", ]
p3 <- ggplot(res, aes(x = Fall_on, y = Functional_score, fill = Fall_on))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("Membrane")+
  ylab("Functional Score")+
  theme_classic2()
p3

res <- fread("../7.eachpro/pSTY.Mitochondrion.position.csv", header = T, stringsAsFactors = F, data.table = F, sep = ",")
res <- res[res$Fall_on != "NaN", ]
p4 <- ggplot(res, aes(x = Fall_on, y = Functional_score, fill = Fall_on))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("Mitochondrion")+
  ylab("Functional Score")+
  theme_classic2()
p4

res <- fread("../7.eachpro/pSTY.Secreted.position.csv", header = T, stringsAsFactors = F, data.table = F, sep = ",")
res <- res[res$Fall_on != "NaN", ]
p5 <- ggplot(res, aes(x = Fall_on, y = Functional_score, fill = Fall_on))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("Secreted")+
  ylab("Functional Score")+
  theme_classic2()
p5

pdf("../8.results/pSTY.function.score.pdf", width = 25, height = 4)
grid.arrange(p1,p2,p3,p4,p5, nrow = 1)
dev.off()

#################################Mutation effects#################################
########mutation in delta regions########
sigregions <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Cytosol.regions", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
sigregions <- sigregions[sigregions$V1 %in% hsdataset$V1, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/mut.point.Cytosol.region.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Membrane.regions", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Membrane", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
sigregions <- sigregions[sigregions$V1 %in% hsdataset$V1, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/mut.point.Membrane.region.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Mitochondrion.regions", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
sigregions <- sigregions[sigregions$V1 %in% hsdataset$V1, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/mut.point.Mitochondrion.region.txt", col.names = T, row.names = F, sep = '\t', quote = F)

sigregions <- fread("../7.eachpro/hs.CPT.delta.scores.regions.sigdelta.Secreted.regions", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Secreted", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
sigregions <- sigregions[sigregions$V1 %in% hsdataset$V1, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../8.results/mut.point.Secreted.region.txt", col.names = T, row.names = F, sep = '\t', quote = F)


res <- fread("../8.results/mut.point.Cytosol.region.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot1 <- melt(res, id.vars = c("ID", "Tumor"))
p1 <- ggplot(res4plot1, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (Cytosol)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/mut.point.Membrane.region.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot1 <- melt(res, id.vars = c("ID", "Tumor"))
p2 <- ggplot(res4plot1, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (Membrane)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/mut.point.Mitochondrion.region.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot1 <- melt(res, id.vars = c("ID", "Tumor"))
p3 <- ggplot(res4plot1, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (Mitochondrion)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

res <- fread("../8.results/mut.point.Secreted.region.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot1 <- melt(res, id.vars = c("ID", "Tumor"))
p4 <- ggplot(res4plot1, aes(x = Tumor, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (Secreted)")+
  theme_linedraw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("../8.results/mut.point.distibute.conbined.boxplot.pdf", width = 10, height = 16)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()

########define best selection########
cptdata1 <- fread("../7.eachpro/hs.CPT.delta.scores.Nucleus", header =T, stringsAsFactors = F, data.table = F)
allcuts <- seq(0,0.9,0.01)

allperc <- c()
for (cuts in allcuts) {
  tmpdat1 <- cptdata1[as.numeric(cptdata1$V8) > cuts, ]
  tmpdat2 <- cptdata1[as.numeric(cptdata1$V8) < -cuts, ]
  allperc <- c(allperc, (nrow(tmpdat1) + nrow(tmpdat2))/nrow(cptdata1))
}

res <- data.frame(cutoff = allcuts,
                  perc = allperc*100)
selres <- res[res$cutoff == 0.1, ]

p <- ggplot(res, aes(x = cutoff, y = perc))+
  geom_line()+
  geom_point(data = selres, aes(x = cutoff, y = perc))+
  geom_text_repel(data = selres, aes(x = cutoff, y = perc), label = "(0.1, 5.17%)")+
  xlab("Delta score")+
  ylab("Significant site percents (%)")+
  theme_classic2()

p
ggsave("../8.results/01.sel.cutoff.pdf", width = 5, height = 4)


allnums <- rnorm(1000)
res <- data.frame(ID = allnums,
                  Num = allnums)

#gg <- res[,list(x=density(res$Num)$x, y=density(res$Num)$y)]

gg <- data.frame(x = density(res$Num)$x,
                 y = density(res$Num)$y)

p<-ggplot(res, aes(x=Num)) +
  geom_density()+
  geom_vline(xintercept = c(-2, 2), lty = "dashed")+
  geom_ribbon(data=subset(gg,x > 2),
              aes(x=x,ymax=y),ymin=0,fill="red", alpha=0.5)+
  geom_ribbon(data=subset(gg,x < -2),
              aes(x=x,ymax=y),ymin=0,fill="blue", alpha=0.5)+
  theme_classic2()

p
ggsave("../8.results/01.sel.cutoff.density.pdf", width = 5, height = 4)



mut_nuc <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
names(mut_nuc) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_nuc$Delta <- mut_nuc$Mut_Score - mut_nuc$Origin_Score
mut_nuc$Class <- "Nucleus"

mut_cyt <- fread("../4.panMut/merged.mut.maf.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)
names(mut_cyt) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_cyt$Delta <- mut_cyt$Mut_Score - mut_cyt$Origin_Score
mut_cyt$Class <- "Cytosol"

mut_mem <- fread("../4.panMut/merged.mut.maf.score.Membrane", header = F, stringsAsFactors = F, data.table = F)
names(mut_mem) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_mem$Delta <- mut_mem$Mut_Score - mut_mem$Origin_Score
mut_mem$Class <- "Membrane"

mut_mit <- fread("../4.panMut/merged.mut.maf.score.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
names(mut_mit) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_mit$Delta <- mut_mit$Mut_Score - mut_mit$Origin_Score
mut_mit$Class <- "Mitochondrion"

mut_sec <- fread("../4.panMut/merged.mut.maf.score.Secreted", header = F, stringsAsFactors = F, data.table = F)
names(mut_sec) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_sec$Delta <- mut_sec$Mut_Score - mut_sec$Origin_Score
mut_sec$Class <- "Secreted"

#mutmerge <- rbind(mut_nuc, mut_cyt, mut_mem, mut_mit, mut_sec)
mutmerge <- mut_nuc
mutmerge$ID <- paste(mutmerge$gene, mutmerge$Amino_Acid_Change, sep = "_")

allcuts <- seq(0,0.9,0.01)
allperc <- c()
for (cuts in allcuts) {
  genesmerge <- mutmerge[abs(mutmerge$Delta) > (cuts/5), ]
  allids <- c(genesmerge[genesmerge$Class == "Nucleus", ]$ID, 
              genesmerge[genesmerge$Class == "Cytosol", ]$ID, 
              genesmerge[genesmerge$Class == "Membrane", ]$ID, 
              genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > (cuts+0.4), ]$ID, 
              genesmerge[genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > (cuts+0.4), ]$ID)
  allperc <- c(allperc, length(unique(allids)) / length(unique(mutmerge$ID)))
}

res <- data.frame(cutoff = allcuts,
                  perc = allperc*100)
selres <- res[res$cutoff == 0.1, ]


p <- ggplot(res, aes(x = cutoff, y = perc))+
  geom_line()+
  geom_point(data = selres, aes(x = cutoff, y = perc))+
  geom_text_repel(data = selres, aes(x = cutoff, y = perc), label = "(0.1, 5.70%)")+
  xlab("Delta score")+
  ylab("Mutation percents (%)")+
  theme_classic2()

p
selres

ggsave("../8.results/01.sel.cutoff.mut.pdf", width = 5, height = 4)



########location sig delta coverage########
mut_nuc <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
names(mut_nuc) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_nuc$Delta <- mut_nuc$Mut_Score - mut_nuc$Origin_Score
mut_nuc$Class <- "Nucleus"

mut_cyt <- fread("../4.panMut/merged.mut.maf.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)
names(mut_cyt) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_cyt$Delta <- mut_cyt$Mut_Score - mut_cyt$Origin_Score
mut_cyt$Class <- "Cytosol"

mut_mem <- fread("../4.panMut/merged.mut.maf.score.Membrane", header = F, stringsAsFactors = F, data.table = F)
names(mut_mem) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_mem$Delta <- mut_mem$Mut_Score - mut_mem$Origin_Score
mut_mem$Class <- "Membrane"

mut_mit <- fread("../4.panMut/merged.mut.maf.score.Mitochondrion", header = F, stringsAsFactors = F, data.table = F)
names(mut_mit) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_mit$Delta <- mut_mit$Mut_Score - mut_mit$Origin_Score
mut_mit$Class <- "Mitochondrion"

mut_sec <- fread("../4.panMut/merged.mut.maf.score.Secreted", header = F, stringsAsFactors = F, data.table = F)
names(mut_sec) <- c("gene", "chrom", "start", "end", "ref", "alt", "Amino_Acid_Change", "effect", "filter", "UniID", "Origin_Score", "Mut_Score")
mut_sec$Delta <- mut_sec$Mut_Score - mut_sec$Origin_Score
mut_sec$Class <- "Secreted"

mutmerge <- rbind(mut_nuc, mut_cyt, mut_mem, mut_mit, mut_sec)
mutmerge$ID <- paste(mutmerge$gene, mutmerge$Amino_Acid_Change, sep = "_")
genesmerge <- mutmerge[abs(mutmerge$Delta) > 0.15, ]
length(unique(genesmerge$ID))
genesmerge$ID <- paste(genesmerge$gene, genesmerge$Amino_Acid_Change, sep = "_")
genesmerge <- genesmerge[, c("ID", "gene", "Origin_Score", "Mut_Score", "Delta", "Class")]
genesmerge <- unique(genesmerge)
write.table(genesmerge, "../4.panMut/sig.delta.merged.txt", row.names = F, col.names = T, quote = F, sep = "\t")

cags <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
enzymes <- fread("/home/yukai/data/Datasets/HPA/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F)
res <- as.data.frame(table(genesmerge$ID))
names(res) <- c("ID", "Freq")
res$Gene <- gsub("_.*", "", res$ID)
res <- res[res$Gene %in% c(cags$Gene, enzymes$Gene), ]

sel_muts <- c("NFE2L1_p.R674C", "PTEN_p.R14M", "PTEN_p.K13E", "RC3H1_p.T9M", "SLC9A1_p.G6D", 
              "CHFR_p.P255L", "SYNPO2_p.G4E", "DAO_p.V4E", "KAT8_p.R144C", "CMAS_p.R201W",
              "DEAF1_p.K304N", "MMP19_p.L6M")

res4plot <- mutmerge[mutmerge$ID %in% sel_muts, ]
res4plot$Class <- factor(res4plot$Class, c("Nucleus", "Cytosol", "Membrane", "Mitochondrion", "Secreted"))
p <- ggplot(res4plot, aes(x = Class, y = ID, fill = Delta))+
  geom_tile(size = 0.5, na.rm = T, color = "gray50")+
  scale_fill_gradient2(low="#B6C9F0", high="#E93B81",mid = "white", midpoint = 0.)+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p
ggsave("../8.results/mut.selpoint.heatmap.pdf", p, width = 4, height = 5.2)

########sig delta upset########
mutnuc <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)

library(UpSetR)
listInput <- list(Nucleus = genesmerge[genesmerge$Class == "Nucleus", ]$ID, 
                  Cytosol = genesmerge[genesmerge$Class == "Cytosol", ]$ID, 
                  Membrane = genesmerge[genesmerge$Class == "Membrane", ]$ID, 
                  Mitochondrion = genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8, ]$ID, 
                  Secreted = genesmerge[genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6, ]$ID)

upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")

pdf("../8.results/mut.point.multilocs.upset.pdf", width = 6, height = 4)
upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")
dev.off()

########sig delta oncoscore########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
#genesmerge <- genesmerge[(genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8) |
#                           (genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6), ]
genesmerge <- genesmerge[abs(genesmerge$Delta) > 0.4, ]
res <- as.data.frame(table(genesmerge$ID))
names(res) <- c("ID", "Freq")
res$Gene <- gsub("_.*", "", res$ID)
res <- res[res$Freq > 1, ]
length(unique(res$Gene))

e3s <- unique(res$Gene)
allfiles <- fread("/home/yukai/softwares/2020plus-master/tumor_sample_stat.txt", header = F, stringsAsFactors = F, data.table = F)
aikeids <- allfiles$V1

finalres <- data.frame()
for (ids in aikeids) {
  scores <- fread(paste("/home/yukai/softwares/2020plus-master/output_res/", ids, "/output/results/r_random_forest_prediction.txt", sep = ""), header = T, stringsAsFactors = 1, data.table = F)
  scores <- scores[, c(1, 26, 27, 37, 38)]
  scores$Class <- ids
  finalres <- rbind(finalres, scores)
}

##uniq and densityplot
res4plot <- finalres
res4plot <- res4plot[res4plot$gene %in% e3s, ]
res4plot1 <- res4plot[order(res4plot$oncogene, decreasing = T), ]
res4plot1 <- res4plot1 %>% distinct(gene, .keep_all = TRUE)
res4plot2 <- res4plot[order(res4plot$tsg, decreasing = T), ]
res4plot2 <- res4plot2 %>% distinct(gene, .keep_all = TRUE)
res4plot3 <- res4plot[order(res4plot$`driver p-value`), ]
res4plot3 <- res4plot3 %>% distinct(gene, .keep_all = TRUE)
res4plot1 <- res4plot1[order(res4plot1$gene), ]
res4plot2 <- res4plot2[order(res4plot2$gene), ]
res4plot3 <- res4plot3[order(res4plot3$gene), ]
res4plot <- data.frame(gene = res4plot1$gene,
                       oncogene = res4plot1$`oncogene score`,
                       tsg = res4plot2$`tsg score`,
                       pvalue = res4plot3$`driver p-value`)

res4plot$group <- ifelse(res4plot$pvalue < 0.05, "Significant", "None")
res4plot4sig <- res4plot[res4plot$pvalue < 0.01, ]
p1 <- ggplot(res4plot, aes(x = tsg, y = oncogene))+
  geom_point(aes(x = tsg, y = oncogene, color = group), size = 0.6)+
  stat_density_2d(colour="black", lwd = 0.4)+
  scale_color_manual(values = c("gray", "red"))+
  xlab("TSG score")+
  ylab("Oncogene score")+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  geom_text_repel(data = res4plot4sig, aes(label = gene), color = "red", size = 2)+
  theme_bw()
p1
ggsave("../8.results/mut.point.selgenes.oncoscores.pdf", p1, width = 4, height = 2.8)

########sig delta CERES########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
#genesmerge <- genesmerge[(genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8) |
#                           (genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6), ]
genesmerge <- genesmerge[abs(genesmerge$Delta) > 0.4, ]
res <- as.data.frame(table(genesmerge$ID))
names(res) <- c("ID", "Freq")
res$Gene <- gsub("_.*", "", res$ID)
res <- res[res$Freq > 1, ]
length(unique(res$Gene))

e3s <- unique(res$Gene)

allfiles <- fread("/home/yukai/softwares/2020plus-master/tumor_sample_stat.txt", header = F, stringsAsFactors = F, data.table = F)
aikeids <- allfiles$V1
finalres <- data.frame()
for (ids in aikeids) {
  scores <- fread(paste("/home/yukai/softwares/2020plus-master/output_res/", ids, "/output/results/r_random_forest_prediction.txt", sep = ""), header = T, stringsAsFactors = 1, data.table = F)
  scores <- scores[, c(1, 26, 27, 37, 38)]
  scores$Class <- ids
  finalres <- rbind(finalres, scores)
}

##uniq and densityplot
res4plot <- finalres
res4plot <- res4plot[res4plot$gene %in% e3s, ]
res4plot1 <- res4plot[order(res4plot$oncogene, decreasing = T), ]
res4plot1 <- res4plot1 %>% distinct(gene, .keep_all = TRUE)
res4plot2 <- res4plot[order(res4plot$tsg, decreasing = T), ]
res4plot2 <- res4plot2 %>% distinct(gene, .keep_all = TRUE)
res4plot3 <- res4plot[order(res4plot$`driver p-value`), ]
res4plot3 <- res4plot3 %>% distinct(gene, .keep_all = TRUE)
res4plot1 <- res4plot1[order(res4plot1$gene), ]
res4plot2 <- res4plot2[order(res4plot2$gene), ]
res4plot3 <- res4plot3[order(res4plot3$gene), ]
res4plot <- data.frame(gene = res4plot1$gene,
                       oncogene = res4plot1$`oncogene score`,
                       tsg = res4plot2$`tsg score`,
                       pvalue = res4plot3$`driver p-value`)

res4plot$group <- ifelse(res4plot$pvalue < 0.01, "Significant", "None")
siggenes <- unique(res4plot[res4plot$group == "Significant", ]$gene)

ceres <- fread("/home/yukai/data/Datasets/DeMap/CCLE_D2_combined_gene_dep_scores.csv", header = T, stringsAsFactors = F, data.table = F)
ceres$V1 <- gsub(" .*", "", ceres$V1)
rownames(ceres) <- ceres$V1

alltiss <- as.data.frame(table(sub("[^_]*_", "", names(ceres))))
alltiss <- alltiss[alltiss$Freq > 20, ]

e3ceres <- ceres[intersect(rownames(ceres), e3s), ]
e3sval <- data.frame(ID = rownames(e3ceres),
                     Median = rowMedians(as.matrix(e3ceres[, -1]), na.rm = T))
e3sval <- e3sval[order(e3sval$Median), ]

tmpceres <- ceres[, -1]
ceres_md <- data.frame(ID = rownames(tmpceres),
                       Score = rowMeans(tmpceres, na.rm = T),
                       Class = "All")
wca <- ceres_md
wca$Class <- ifelse(wca$ID %in% e3sval$ID[1:30], "0.Sigpros", "1.Others")
wca <- wca[is.nan(wca$Score) == FALSE, ]
library(plyr)
mu <- ddply(wca, "Class", summarise, grp.mean=mean(Score))
p<-ggplot(wca, aes(x=Score, fill=Class)) +
  geom_density(alpha=0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Class),
             linetype="dashed")+
  theme_classic2()
print(p)
ggsave("../8.results/mut.point.selgenes.ceres.pdf", p, width = 6, height = 5)


########phosmimic########
nlsmimic <- fread("../4.panMut/hs.CPT.delta.scores.regions.sigdelta.Nucleus.regions.mutAA.results", 
                  header = T, stringsAsFactors = F, data.table = F)
names(nlsmimic) <- c("UniID", "Gene", "Seq", "Score", "Region", "toA", "toD", "ScoretoA", "ScoretoD")

res <- nlsmimic[, c(1, 2, 4, 9)]
res4plot <- melt(res, id.vars = c("UniID","Gene"))

p1 <- ggplot(res4plot, aes(x = variable, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means()+
  xlab("")+
  ylab("Predicted Score")+
  theme_classic2()
p1

ggsave("../8.results/mut.point.phos.nlstoD.pdf", p1, width = 4, height = 3)

cor.test(nlsmimic$ScoretoA - nlsmimic$Score, nlsmimic$ScoretoD - nlsmimic$Score)

nrow(nlsmimic[nlsmimic$Score == nlsmimic$ScoretoD, ])

#summary
nlsmic <- nlsmimic
nlsmic$Class <- "NoChange"
nlsmic$Class <- ifelse(nlsmic$ScoretoD - nlsmic$Score > 0.05, "Increased", nlsmic$Class)
nlsmic$Class <- ifelse(nlsmic$ScoretoD - nlsmic$Score < -0.05, "Decreased", nlsmic$Class)
nlsmic$Class <- ifelse(nlsmimic$Score == nlsmimic$ScoretoD, "NoSites", nlsmic$Class)
table(nlsmic$Class)
res <- as.data.frame(table(nlsmic$Class))
res$Perc <- round((res$Freq / sum(res$Freq))*100, 2)

df2 <- res %>% 
  mutate(csum = rev(cumsum(rev(Freq))), 
         pos = Freq/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Freq/2, pos))

p <- ggplot(res, aes(x = "" , y = Freq, fill = Var1)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(Perc, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()

ggsave("../8.results/mut.point.phos.nlstoD.distrib.pdf", p, width = 4, height = 3)

kinsub <- readgmt(readLines("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/enzyme.substrate.psty.gmt"))
kinsub$Subname <- gsub("#.*", "", kinsub$genes)
kinsub <- unique(kinsub[, c(1,3)])
sel_mimic <- nlsmimic[(nlsmimic$ScoretoD > nlsmimic$Score) & (nlsmimic$ScoretoA - nlsmimic$Score < -0.1), ]
res_hall <- enricher(sel_mimic$Gene, pvalueCutoff = 1,
                     TERM2GENE = kinsub, qvalueCutoff = 1, minGSSize = 1)
res2subtop15 <- res_hall@result
res2subtop15$link <- -log10(res2subtop15$pvalue)
res2subtop15 <- res2subtop15[res2subtop15$pvalue < 0.05, ]
res2subtop15$ID <- factor(res2subtop15$ID, rev(res2subtop15$ID))
p2 <- ggplot(res2subtop15, aes(x = ID, y = link))+
  geom_bar(stat = "identity", fill = "lightblue", color = "black", width = 0.7)+
  theme_classic2()+
  #scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()
p2
ggsave("../8.results/mut.nls.psty.upkinase.pdf", p2, width = 5, height = 3)

##all kins survival
kinsub <- readgmt(readLines("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/enzyme.substrate.psty.gmt"))
kinsub$Subname <- gsub("#.*", "", kinsub$genes)
kinsub <- unique(kinsub[, c(1,3)])
sel_mimic <- nlsmimic[(nlsmimic$ScoretoD > nlsmimic$Score) & (nlsmimic$ScoretoA - nlsmimic$Score < -0.1), ]
res_hall <- enricher(sel_mimic$Gene, pvalueCutoff = 1,
                     TERM2GENE = kinsub, qvalueCutoff = 1, minGSSize = 1)
res2subtop15 <- res_hall@result
res2subtop15 <- res2subtop15[res2subtop15$pvalue < 0.05, ]
allkins <- res2subtop15$ID

tumos <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)

finalres <- data.frame()
for (selkins in allkins) {
  for (selcancer in tumos$V1) {
    print(paste(selkins, "in", selcancer, "started!!!!", sep = " "))
    expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", selcancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
    rownames(expmtr) <- expmtr$Ensembl_ID
    expmtr <- expmtr[, -1]
    expmtr <- expmtr[, as.numeric(substr(names(expmtr), 14, 15)) < 10]
    survinfo <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", selcancer, ".survival.tsv", sep = ""), 
                      header = T, stringsAsFactors = F, data.table = F)
    rownames(survinfo) <- survinfo$sample
    tmpres <- surv_bestcut(expmtr, selkins, survinfo, num = 20)
    tmpres$Class <- selcancer
    finalres <- rbind(finalres, tmpres[order(tmpres$Pvalue_OS), ][1, ])
  }
}

write.table(finalres, "../8.results/mut.point.phos.kinases.survival.txt", row.names = F, col.names = T, quote = F, sep = "\t")
finalres <- fread("../8.results/mut.point.phos.kinases.survival.txt", header = T, stringsAsFactors = F, data.table = F)
finalres$HR_OS <- log(finalres$HR_OS)

res4sig <- finalres[finalres$Pvalue_OS < 0.05, ]
finalres$ID <- factor(finalres$ID, rev(allkins))
p <- ggplot()+
  geom_tile(data = finalres, aes(x = Class, y = ID, fill = HR_OS))+
  geom_point(data = res4sig, aes(x = Class, y = ID), shape= 8)+
  theme_classic2()+
  scale_fill_gradient2(low="blue", high="red",mid = "white", midpoint = 0, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p
ggsave("../8.results/mut.point.phos.kinases.survival.pdf", width = 5, height = 3)




########selected gene survival########
tumos <- fread("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/tumor_sample_stat.txt", header = F, stringsAsFactors = F, data.table = F)
selgene <- "PLK3"

for (selcancer in tumos$V1) {
  expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", selcancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                  header = T, stringsAsFactors = F, data.table = F)
  rownames(expmtr) <- expmtr$Ensembl_ID
  expmtr <- expmtr[, -1]
  keep <- rowSums(expmtr > 1) >= ncol(expmtr)/2 #a Count>0 in at least 3 samples
  expmtr <- expmtr[keep,]
  expmtr <- expmtr[, as.numeric(substr(names(expmtr), 14, 15)) < 10]
  survinfo <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", selcancer, ".survival.tsv", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
  rownames(survinfo) <- survinfo$sample
  commonsams <- intersect(rownames(survinfo), names(expmtr))
  survinfo <- survinfo[commonsams, ]
  expmtr <- expmtr[, commonsams]
  cluster_surv <- survinfo
  sevalue <- as.numeric(expmtr[selgene, ])
  cluster_surv$Type = ifelse(sevalue > median(sevalue), "0.High", "1.Low")
  tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
  fit <- survfit(Surv(OS.time, OS) ~ Type, data = cluster_surv)
  os <- survplot(cluster_surv, type = paste("OS", selcancer, sep = "-"), fit = fit, pval = tmp$logtest[3])
  print(os$plot)
}

selgene <- "PLK3"
selcancer <- "LGG"
expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", selcancer, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]
keep <- rowSums(expmtr > 1) >= ncol(expmtr)/2 #a Count>0 in at least 3 samples
expmtr <- expmtr[keep,]
expmtr <- expmtr[, as.numeric(substr(names(expmtr), 14, 15)) < 10]
survinfo <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", selcancer, ".survival.tsv", sep = ""), 
                  header = T, stringsAsFactors = F, data.table = F)
rownames(survinfo) <- survinfo$sample
commonsams <- intersect(rownames(survinfo), names(expmtr))
survinfo <- survinfo[commonsams, ]
expmtr <- expmtr[, commonsams]
cluster_surv <- survinfo
sevalue <- as.numeric(expmtr[selgene, ])
cluster_surv$Type = ifelse(sevalue > median(sevalue), "0.High", "1.Low")
tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
fit <- survfit(Surv(OS.time, OS) ~ Type, data = cluster_surv)
os <- survplot(cluster_surv, type = paste("OS", selcancer, sep = "-"), fit = fit, pval = tmp$logtest[3])
print(os$plot)

pdf("../8.results/mut.nls.upkinase.surv.pdf", width = 6, height = 3.5)
print(os$plot)
dev.off()

#################################Mut Analysis#################################
########TCGA mutation count########
alltumors <- fread("/home/yukai/data/Datasets/TCGA20201022/TCGA.samples.txt",
                   header = F, stringsAsFactors = F, data.table = F)
allcounts <- fread("/home/yukai/data/Datasets/TCGA20201022/mutation/mutect2/tmp.txt",
                   header = F, stringsAsFactors = F, data.table = F)
res <- cbind(alltumors, allcounts)
names(res) <- c("Tumor", "Count")

#50000以上突变

########site&gene&enrichment analysis########
mutnuc <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)

library(UpSetR)
listInput <- list(Nucleus = genesmerge[genesmerge$Class == "Nucleus", ]$ID, 
                  Cytosol = genesmerge[genesmerge$Class == "Cytosol", ]$ID, 
                  Membrane = genesmerge[genesmerge$Class == "Membrane", ]$ID, 
                  Mitochondrion = genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8, ]$ID, 
                  Secreted = genesmerge[genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6, ]$ID)

upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")

pdf("../8.results/mut.point.multilocs.upset.pdf", width = 6, height = 4)
upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")
dev.off()

listInput <- list(Nucleus = genesmerge[genesmerge$Class == "Nucleus", ]$gene, 
                  Cytosol = genesmerge[genesmerge$Class == "Cytosol", ]$gene, 
                  Membrane = genesmerge[genesmerge$Class == "Membrane", ]$gene, 
                  Mitochondrion = genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8, ]$gene, 
                  Secreted = genesmerge[genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6, ]$gene)

upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")

pdf("../8.results/mut.point.multilocs.gene.upset.pdf", width = 6, height = 4)
upset(fromList(listInput), order.by = "freq", 
      sets.bar.color = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF","#B05E27FF", "#D4AC2BFF"),
      matrix.color = "#4DBBD5FF")
dev.off()

#####enriched pathways
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

ewp1 <- enricher(genesmerge[genesmerge$Class == "Nucleus", ]$gene, pvalueCutoff = 1,
                     TERM2GENE = gobp, qvalueCutoff = 1, minGSSize = 1)
ewp2 <- enricher(genesmerge[genesmerge$Class == "Cytosol", ]$gene, pvalueCutoff = 1,
                 TERM2GENE = gobp, qvalueCutoff = 1, minGSSize = 1)
ewp3 <- enricher(genesmerge[genesmerge$Class == "Membrane", ]$gene, pvalueCutoff = 1,
                 TERM2GENE = gobp, qvalueCutoff = 1, minGSSize = 1)
ewp4 <- enricher(genesmerge[genesmerge$Class == "Mitochondrion", ]$gene, pvalueCutoff = 1,
                 TERM2GENE = gobp, qvalueCutoff = 1, minGSSize = 1)
ewp5 <- enricher(genesmerge[genesmerge$Class == "Secreted", ]$gene, pvalueCutoff = 1,
                 TERM2GENE = gobp, qvalueCutoff = 1, minGSSize = 1)

path1 <- ewp1@result[ewp1@result$pvalue < 0.01, ][1:10, ]
path2 <- ewp2@result[ewp2@result$pvalue < 0.01, ][1:10, ]
path3 <- ewp3@result[ewp3@result$pvalue < 0.01, ][1:10, ]
path4 <- ewp4@result[ewp4@result$pvalue < 0.01, ][1:10, ]
path5 <- ewp5@result[ewp5@result$pvalue < 0.01, ][1:10, ]
pathes <- Reduce(union,list(
  "Nucleus" = path1$ID,
  "Cytosol" = path2$ID,
  "Membrane" = path3$ID,
  "Mitochondrion" = path4$ID,
  "Secreted" = path5$ID))

new_path1 <- ewp1@result[pathes, ]
new_path2 <- ewp2@result[pathes, ]
new_path3 <- ewp3@result[pathes, ]
new_path4 <- ewp4@result[pathes, ]
new_path5 <- ewp5@result[pathes, ]
upper_tri <- data.frame(ID = pathes,
                        Nucleus = -log10(new_path1$pvalue),
                        Cytosol = -log10(new_path2$pvalue),
                        Membrane = -log10(new_path3$pvalue),
                        Mitochondrion = -log10(new_path4$pvalue),
                        Secreted = -log10(new_path5$pvalue))

upper_tri$Nucleus <- ifelse(upper_tri$ID %in% path1$ID, upper_tri$Nucleus, 0)
upper_tri$Cytosol <- ifelse(upper_tri$ID %in% path2$ID, upper_tri$Cytosol, 0)
upper_tri$Membrane <- ifelse(upper_tri$ID %in% path3$ID, upper_tri$Membrane, 0)
upper_tri$Mitochondrion <- ifelse(upper_tri$ID %in% path4$ID, upper_tri$Mitochondrion, 0)
upper_tri$Secreted <- ifelse(upper_tri$ID %in% path5$ID, upper_tri$Secreted, 0)

library(forcats)
upper_tri_melt <- reshape2::melt(upper_tri, na.rm = TRUE, id.vars = 'ID')
names(upper_tri_melt) <- c('Pathway', 'SubTypes', 'Values')
upper_tri_melt$Values <- log10(upper_tri_melt$Values+1)
upper_tri_melt$Pathway <- fct_inorder(upper_tri_melt$Pathway)

pdf('../8.results/mut.analy.deltagene.pathways.pdf', width = 15, height = 12)
ggplot(data=upper_tri_melt, aes(SubTypes, y=Pathway, fill=Values))+
  geom_tile(color="black", size = 0.5)+
  scale_fill_gradient2(low = "white", high = "red", space = "Lab", name="-logFDR")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_flip()
dev.off()

########nuclear top proteins########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
library(UpSetR)
deltapros <- unique(genesmerge[genesmerge$Class == "Nucleus" & abs(genesmerge$Delta) > 0.2, ]$gene)
deltapros <- deltapros[grep("-", deltapros, invert = T)]

hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
ewp1 <- enricher(deltapros, pvalueCutoff = 1,
                 TERM2GENE = hall, qvalueCutoff = 1, minGSSize = 1)
View(ewp1@result)

res <- ewp1@result[, c(1, 8)]
res <- res %>% separate_rows(geneID, sep = "/")
res <- res %>% distinct(geneID, .keep_all = TRUE)
res <- res[res$ID %in% ewp1@result$ID[c(1:9, 12, 13)], ]

enzymes <- fread("/home/yukai/data/Datasets/HPA/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F)
transcr <- fread("/home/yukai/data/Datasets/HPA/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
fdas <- fread("/home/yukai/data/Datasets/HPA/protein_class_FDA.tsv", header = T, stringsAsFactors = F, data.table = F)
cags <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
metastasis <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/EMT.Meta.Metabolic.oncoTSG.list", header = F, stringsAsFactors = F, data.table = F)
nucrec <- fread("/home/yukai/data/Datasets/HPA/protein_class_NuclearReceptors.tsv", header = T, stringsAsFactors = F, data.table = F)
knownnuc <- fread("/home/yukai/projects/CBAmodel/subCellLoc/0.data/0.2Nucleus_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
uni2pro <- fread("/home/yukai/projects/CBAmodel/subCellLoc/0.dataset/uniprot-human.tab", header = T, stringsAsFactors = F, data.table = F)
uni2pro$`Gene names` <- gsub(" .*", "", uni2pro$`Gene names`)
knownnuc <- merge(knownnuc, uni2pro, by.x = "V1", by.y = "Entry")

res$Enzymes <- ifelse(res$geneID %in% enzymes$Gene, "Yes", "No")
res$TFs <- ifelse(res$geneID %in% transcr$Gene, "Yes", "No")
res$FDAs <- ifelse(res$geneID %in% fdas$Gene, "Yes", "No")
res$CAGs <- ifelse(res$geneID %in% cags$Gene, "Yes", "No")
res$Metastasis <- ifelse(res$geneID %in% metastasis[metastasis$V2 == "Meta", ]$V3, "Yes", "No")
res$NuclearRece <- ifelse(res$geneID %in% nucrec$Gene, "Yes", "No")
res$NucLoc <- ifelse(res$geneID %in% knownnuc$`Gene names`, "Yes", "No")

write.table(res, "../8.results/mut.analy.nucleus.genes.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cancerhalls <- fread("../0.dataset/cancerhallmarks.txt", header = T, stringsAsFactors = F, data.table = F)
cancerhalls <- cancerhalls[cancerhalls$Symbol %in% res$geneID, ]

res4anno <- melt(res, id.vars = c("ID", "geneID"))
res4anno <- res4anno[res4anno$value == "Yes" & res4anno$variable != "NucLoc", ]
write.table(res4anno, "../8.results/mut.analy.nucleus.genes.info.txt", row.names = F, col.names = T, sep = "\t", quote = F)


########nuclear gene and survival########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
deltapros <- unique(genesmerge[genesmerge$Class == "Secreted", ])
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
selnum <- 10
mutsamnum = c()
othnum = c()
mutVsctrSurv = c()
upper <- c()
lower <- c()
mutVsctrSurvPva = c()
tumortype = c()
for (cancer in tumorlist$V1) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut$ID <- paste(singmut$gene, singmut$Amino_Acid_Change, sep = "_")
  survfile <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", cancer, ".survival.tsv", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
  survfile$OS.days <- survfile$OS.time
  survfile$OS.status <- survfile$OS
  
  singmut <- singmut[singmut$ID %in% deltapros$ID, ]
  singmut4res <- as.data.frame(table(singmut$Sample_ID))
  mutsam <- as.character(singmut4res[singmut4res$Freq > selnum, ]$Var1)
  othsam <- setdiff(unique(singmut$Sample_ID), mutsam)
  othnum <- c(othnum, length(othsam))
  if (length(mutsam) < 3) {
    mutsamnum = c(mutsamnum, length(mutsam))
    mutVsctrSurv = c(mutVsctrSurv, NA)
    mutVsctrSurvPva = c(mutVsctrSurvPva, NA)
    upper <- c(upper, NA)
    lower <- c(lower, NA)
    tumortype = c(tumortype, cancer)
  }else{
    mutsamnum = c(mutsamnum, length(mutsam))
    tumortype = c(tumortype, cancer)
    survfile$Class <- NA
    survfile$Class = ifelse(survfile$sample %in% mutsam, "0.mut", survfile$Class)
    survfile$Class = ifelse(survfile$sample %in% othsam, "1.non", survfile$Class)
    survfile <- na.omit(survfile)
    if (nrow(survfile[survfile$Class == "0.mut", ]) < 2) {
      mutVsctrSurv = c(mutVsctrSurv, 1)
      mutVsctrSurvPva <- c(mutVsctrSurvPva, 1)
      upper <- c(upper, 1)
      lower <- c(lower, 1)
    }else{
      tmp <- summary(coxph((Surv(OS.days, OS.status)) ~ Class, data = survfile))
      mutVsctrSurv = c(mutVsctrSurv, tmp$conf.int[[1]])
      mutVsctrSurvPva <- c(mutVsctrSurvPva, tmp$logtest[[3]])
      upper <- c(upper, tmp$conf.int[[4]])
      lower <- c(lower, tmp$conf.int[[3]])
      fit <- survfit(Surv(OS.days, OS.status) ~ Class, data = survfile)
      os <- survplot(survfile, type = paste("OS", cancer, sep = "-"), fit = fit, pval = tmp$logtest[[3]])
      print(os)
      pdf(paste("../8.results/mut.analy.nucleus.", cancer, ".survival.pdf", sep = ""), width = 5, height = 3)
      print(os$plot)
      dev.off()
    }
  }
}

res <- data.frame(Num = mutsamnum,
                  OthNum = othnum,
                  Surv = mutVsctrSurv,
                  SurvPva = mutVsctrSurvPva,
                  Upper = upper,
                  Lower = lower,
                  Cancer = tumortype,
                  stringsAsFactors = F)

write.table(res, "../8.results/mut.analy.secreted.surv.txt", row.names = F, col.names = T, quote = F, sep = "\t")

res1 <- fread("../8.results/mut.analy.nucleus.surv.txt", header = T, stringsAsFactors = F, data.table = F)
res1$Class <- "Nucleus"
res2 <- fread("../8.results/mut.analy.cytosol.surv.txt", header = T, stringsAsFactors = F, data.table = F)
res2$Class <- "Cytosol"
res3 <- fread("../8.results/mut.analy.membrane.surv.txt", header = T, stringsAsFactors = F, data.table = F)
res3$Class <- "Membrane"
res4 <- fread("../8.results/mut.analy.mito.surv.txt", header = T, stringsAsFactors = F, data.table = F)
res4$Class <- "Mitochondrion"
res5 <- fread("../8.results/mut.analy.secreted.surv.txt", header = T, stringsAsFactors = F, data.table = F)
res5$Class <- "Secreted"
res <- rbind(res1, res2, res3, res4, res5)
res$logp <- -log10(res$SurvPva)
res4sig <- res[res$SurvPva < 0.1, ]

p <- ggplot()+
  geom_tile(data = res, aes(x = Cancer, y = Class, fill = logp))+
  geom_point(data = res4sig, aes(x = Cancer, y = Class), shape= 8)+
  theme_classic2()+
  scale_fill_gradient2(low="#B6C9F0", high="#E93B81",mid = "#FFE5E2", midpoint = 1, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("../8.results/mut.analy.pancancer.surv.pdf", width = 5, height = 2)

########nuclear gene and immune########
library(GSVA)
library(msigdbr)
library(GSEABase)
tcgaMtr <- fread("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-UCEC.htseq_fpkm.tsv.cv.txt", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
tcgaMtr <- tcgaMtr[, -1]
colData <- data.frame(Type = ifelse(substr(colnames(tcgaMtr), 14, 15) < 10, 'Tumor', 'Normal'),
                      ID = colnames(tcgaMtr))
colData$ID <- as.character(colData$ID)
colData <- colData[colData$Type == 'Tumor', ]
tcgaMtr <- tcgaMtr[, colData$ID]
neoantigen <- read.delim("/home/yukai/data/Datasets/TCGA_old/AnalysisDataset/TCGA_PCA.mc3.v0.2.8.CONTROLLED.filtered.sample_neoantigens_10062017.tsv", header = T, row.names = 1, stringsAsFactors = F)
rownames(neoantigen) <- paste(rownames(neoantigen), "-01A", sep = "")
common_samples <- intersect(colData$ID, rownames(neoantigen))
immuneCirc <- getGmt("/home/yukai/work/gc_data/gc_combine_analysis/20200707_datas/0.database/immune.circle.signature.gmt")
tcgaMtr <- as.matrix(tcgaMtr)
res <- gsva(tcgaMtr, immuneCirc, method = "ssgsea")
res.norm <- res
res.norm[c(1:3), ] <- -res.norm[c(1:3), ]
res.norm <- as.data.frame(res.norm)
res.norm <- res.norm[, common_samples]
neoantigen <- neoantigen[common_samples, ]
res.norm <- rbind(res.norm, neoantigen$neoantigen_num)
colnames(res.norm) = common_samples
rownames(res.norm)[nrow(res.norm)] <- "2.neoantigen"
res.norm <- res.norm[order(rownames(res.norm),decreasing = F), ]

nucpros <- fread("../8.results/mut.analy.nucleus.genes.info.txt", header = T, stringsAsFactors = F, data.table = F)
cancerhalls <- fread("../0.dataset/cancerhallmarks.txt", header = T, stringsAsFactors = F, data.table = F)
cancerhalls <- cancerhalls[cancerhalls$Symbol %in% nucpros$geneID, ]
nucexp <- tcgaMtr[unique(cancerhalls$Symbol), common_samples]

rescor <- two_matrix_cor(as.matrix(t(nucexp)), as.matrix(t(res.norm)))
rescor$logp <- -log10(rescor$pva)
rescor4sig <- rescor[rescor$pva < 0.01, ]

p <- ggplot()+
  geom_tile(data = rescor, aes(x = ID1, y = ID2, fill = logp))+
  geom_point(data = rescor4sig, aes(x = ID1, y = ID2), shape= 8)+
  theme_classic2()+
  scale_fill_gradient2(low="#B6C9F0", high="#E93B81",mid = "#FFE5E2", midpoint = 1, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
ggsave("../8.results/mut.analy.nucleus.ucec.immune.pdf", width = 6, height = 3)



########mito top proteins########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
library(UpSetR)
deltapros <- unique(genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.7, ]$gene)
deltapros <- deltapros[grep("-", deltapros, invert = T)]

hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
ewp1 <- enricher(deltapros, pvalueCutoff = 1,
                 TERM2GENE = hall, qvalueCutoff = 1, minGSSize = 1)
View(ewp1@result)

res <- ewp1@result[, c(1, 8)]
res <- res %>% separate_rows(geneID, sep = "/")
res <- res %>% distinct(geneID, .keep_all = TRUE)
res <- res[res$ID %in% ewp1@result$ID[c(1:6, 8, 10, 11)], ]

enzymes <- fread("/home/yukai/data/Datasets/HPA/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F)
transcr <- fread("/home/yukai/data/Datasets/HPA/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
fdas <- fread("/home/yukai/data/Datasets/HPA/protein_class_FDA.tsv", header = T, stringsAsFactors = F, data.table = F)
cags <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
metastasis <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/EMT.Meta.Metabolic.oncoTSG.list", header = F, stringsAsFactors = F, data.table = F)
nucrec <- fread("/home/yukai/data/Datasets/HPA/protein_class_NuclearReceptors.tsv", header = T, stringsAsFactors = F, data.table = F)
knownnuc <- fread("/home/yukai/projects/CBAmodel/subCellLoc/0.data/0.2Mitochondrion_positive_data.txt", header = F, stringsAsFactors = F, data.table = F)
uni2pro <- fread("/home/yukai/projects/CBAmodel/subCellLoc/0.dataset/uniprot-human.tab", header = T, stringsAsFactors = F, data.table = F)
uni2pro$`Gene names` <- gsub(" .*", "", uni2pro$`Gene names`)
knownnuc <- merge(knownnuc, uni2pro, by.x = "V1", by.y = "Entry")

res$Enzymes <- ifelse(res$geneID %in% enzymes$Gene, "Yes", "No")
res$TFs <- ifelse(res$geneID %in% transcr$Gene, "Yes", "No")
res$FDAs <- ifelse(res$geneID %in% fdas$Gene, "Yes", "No")
res$CAGs <- ifelse(res$geneID %in% cags$Gene, "Yes", "No")
res$Metastasis <- ifelse(res$geneID %in% metastasis[metastasis$V2 == "Meta", ]$V3, "Yes", "No")
res$NuclearRece <- ifelse(res$geneID %in% nucrec$Gene, "Yes", "No")
res$NucLoc <- ifelse(res$geneID %in% knownnuc$`Gene names`, "Yes", "No")

write.table(res, "../8.results/mut.analy.mito.genes.txt", row.names = F, col.names = T, sep = "\t", quote = F)

cancerhalls <- fread("../0.dataset/cancerhallmarks.txt", header = T, stringsAsFactors = F, data.table = F)
cancerhalls <- cancerhalls[cancerhalls$Symbol %in% res$geneID, ]

res4anno <- melt(res, id.vars = c("ID", "geneID"))
res4anno <- res4anno[res4anno$value == "Yes" & res4anno$variable != "NucLoc", ]
write.table(res4anno, "../8.results/mut.analy.mito.genes.info.txt", row.names = F, col.names = T, sep = "\t", quote = F)




#################################PTEN selected#################################
########PTEN NLPT########
cptdata_wt <- fread("../4.panMut/mutplot/wt.res.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
#cptdata_mut <- fread("../4.panMut/mutplot/mut.res.CPT.delta.scores.regions", header =F, stringsAsFactors = F, data.table = F)
cptdata_wt_ori <- fread("../4.panMut/mutplot/wt.res.CPT.delta.scores", header =F, stringsAsFactors = F, data.table = F)
cptdata_mut_ori <- fread("../4.panMut/mutplot/mut.res.CPT.delta.scores", header =F, stringsAsFactors = F, data.table = F)
names(cptdata_wt) <- c("UniID", "Pos", "Score", "Strunc", "Delta")
names(cptdata_mut) <- c("UniID", "Pos", "Score", "Strunc", "Delta")
names(cptdata_wt_ori) <- c("UniID", "Pos", "Score", "Strunc", "Delta")
names(cptdata_mut_ori) <- c("UniID", "Pos", "Score", "Strunc", "Delta")
cptdata_wt$Type <- "WT"
cptdata_mut_ori$Type <- "MUT"
select_cptdata <- cptdata_wt
select_cptdata$Strunc <- cptdata_wt_ori$Strunc
select_cptdata <- rbind(select_cptdata, cptdata_mut_ori)
#select_cptdata$Strunc_mut <- cptdata_mut_ori$Strunc
##single protein version
library(scales)
p2 <- ggplot() +
  geom_line(data = select_cptdata,aes(x = Pos,y = Strunc, linetype = Type), color = "red") +
  geom_line(data = select_cptdata,aes(x = Pos,y = rescale(Delta, c(0, 1)), linetype = Type),
            color = "blue") +
  scale_y_continuous(limits = c(0,1),breaks = c(seq(0,1,0.2)),
                     sec.axis = sec_axis( ~rescale(., c(min(select_cptdata$Delta),
                                                        max(select_cptdata$Delta))),
                                          name = "Delta Coding Probability")) +
  labs(title="Coding Probability Trajectory of protein ")+
  xlab("Protein sequence length")+
  ylab("Coding Probability")+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p2
ggsave("../8.results/mut.before.after.nlpt.pdf", p2, width = 5.5, height = 4)

########others########





########PTEN affected genes########
protein_de_via_wilcox = function(mat) {
  test = mat[names(mat) %in% tmp_clininfo[tmp_clininfo$types == 'Mut', 1]]
  cont = mat[names(mat) %in% tmp_clininfo[tmp_clininfo$types == 'Non', 1]]
  tmp = wilcox.test(as.numeric(test), as.numeric(cont))
  return(tmp$p.value)
}
mutpairs <- c("p.R14M", "p.R14S","p.K13E","p.R11K","p.R15I")

allmuts <- fread("/home/yukai/data/TCGA/mutect_tcga/TCGA-UCEC.mutect2_snv.tsv", header = T, stringsAsFactors = F, data.table = F)
mutsams <- unique(allmuts[allmuts$Amino_Acid_Change %in% mutpairs, ]$Sample_ID)
nonsams <- unique(allmuts[!(allmuts$Amino_Acid_Change %in% mutpairs), ]$Sample_ID)

expmtr <- fread("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-UCEC.htseq_fpkm.tsv.cv.txt", 
                header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]

exp_tes <- expmtr[, names(expmtr) %in% mutsams]
exp_ctr <- expmtr[, names(expmtr) %in% nonsams]

tmp_clininfo <- data.frame(ID = c(mutsams, nonsams),
                           types = c(rep("Mut", length(mutsams)), rep("Non", length(nonsams))))

deres <- data.frame(ID = rownames(expmtr),
                    Meantest = rowMedians(as.matrix(exp_tes)),
                    Meanctr = rowMedians(as.matrix(exp_ctr)),
                    DEFC = log2((rowMedians(as.matrix(exp_tes)) + 1) / (rowMedians(as.matrix(exp_ctr)) + 1)),
                    DEPva = apply(expmtr, 1, protein_de_via_wilcox))


promtr <- fread("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/RPPA/UCEC.pro.matrix.txt", 
                header = T, stringsAsFactors = F, data.table = F)
promtr <- na.omit(promtr)
rownames(promtr) <- promtr$Gene
names(promtr) <- substr(names(promtr), 1, 16)
promtr <- promtr[, -1]
pro_tes <- promtr[, names(promtr) %in% mutsams]
pro_ctr <- promtr[, names(promtr) %in% nonsams]
tmp_clininfo <- data.frame(ID = c(mutsams, nonsams),
                           types = c(rep("Mut", length(mutsams)), rep("Non", length(nonsams))))
depros <- data.frame(ID = rownames(promtr),
                    Meantest = rowMedians(as.matrix(pro_tes)),
                    Meanctr = rowMedians(as.matrix(pro_ctr)),
                    DEFC = rowMedians(as.matrix(pro_tes)) - rowMedians(as.matrix(pro_ctr)),
                    DEPva = apply(promtr, 1, protein_de_via_wilcox))


selgene = "CCND1"

res4plot <- data.frame(ID = c(names(expmtr), names(promtr), names(promtr)),
                       Pros = c(rep("CCND1", ncol(expmtr)), rep("P53", ncol(promtr)), rep("P27_pT198", ncol(promtr))),
                       Class = c(rep("RNA", ncol(expmtr)), rep("protein", ncol(promtr)), rep("protein", ncol(promtr))),
                       ExpValue = c(as.numeric(expmtr["CCND1", ]), as.numeric(promtr["P53", ]), as.numeric(promtr["P27_pT198", ]))
                       )

res4plot <- merge(res4plot, tmp_clininfo, by.x = "ID", by.y = "ID")

p1 <- ggplot(res4plot, aes(x = types, y = ExpValue, fill = types))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", method = "wilcox")+
  ylab("mRNA (protein) expression level")+
  theme_classic2()
q <- facet(p1, facet.by = "Pros", scales = "free_y")
q
pdf("../8.results/01.pten.dn.target.pdf", width = 7, height =4)
print(q)
dev.off()

########select TF and its subs########
protein_de_via_wilcox = function(mat) {
  test = mat[names(mat) %in% tmp_clininfo[tmp_clininfo$types == 'Mut', 1]]
  cont = mat[names(mat) %in% tmp_clininfo[tmp_clininfo$types == 'Non', 1]]
  tmp = wilcox.test(as.numeric(test), as.numeric(cont))
  return(tmp$p.value)
}
tsggenes <- fread("../0.dataset/Human_TSGs.txt", header = T, stringsAsFactors = F, data.table = F, sep = "\t")

nucmut <- fread("../4.panMut/singmut.delta.sig.annote.Nucleus", 
                header = T, stringsAsFactors = F, data.table = F)
nucmut <- nucmut[nucmut$gene %in% tsggenes$GeneSymbol, ]
nucmut <- nucmut[nucmut$deltaScore < -0.1, ]
alltfs <- fread("../0.dataset/TF-Target-information.txt", header = T, stringsAsFactors = F, data.table = F)
nucmut <- nucmut[nucmut$gene %in% alltfs$TF, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)

finalres <- data.frame()
for (tfs in unique(nucmut$gene)) {
  for (tumors in tumorlist$V1) {
    allmuts <- fread(paste("/home/yukai/data/TCGA/mutect_tcga/TCGA-", tumors, ".mutect2_snv.tsv", sep = ""), 
                     header = T, stringsAsFactors = F, data.table = F)
    mutpairs <- nucmut[nucmut$gene == tfs, ]$Amino_Acid_Change
    subgenes <- unique(alltfs[alltfs$TF == tfs, ]$target)
    mutsams <- unique(allmuts[allmuts$Amino_Acid_Change %in% mutpairs, ]$Sample_ID)
    nonsams <- unique(allmuts[!(allmuts$Amino_Acid_Change %in% mutpairs), ]$Sample_ID)
    
    expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", tumors, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
    rownames(expmtr) <- expmtr$Ensembl_ID
    expmtr <- expmtr[, -1]
    expmtr <- expmtr[rownames(expmtr) %in% subgenes, ]
    
    exp_tes <- expmtr[, names(expmtr) %in% mutsams]
    exp_ctr <- expmtr[, names(expmtr) %in% nonsams]
    
    if (length(mutsams) < 5) {
      print("skipped one TF!!!!!!");
    }else{
      print("Useful TF!!!!!!");
      tmp_clininfo <- data.frame(ID = c(mutsams, nonsams),
                                 types = c(rep("Mut", length(mutsams)), rep("Non", length(nonsams))))
      
      deres <- data.frame(ID = rownames(expmtr),
                          TFs = tfs,
                          Cancers = tumors,
                          Meantest = rowMedians(as.matrix(exp_tes)),
                          Meanctr = rowMedians(as.matrix(exp_ctr)),
                          nTest = ncol(exp_tes),
                          nCtr = ncol(exp_ctr),
                          DEFC = log2((rowMedians(as.matrix(exp_tes)) + 1) / (rowMedians(as.matrix(exp_ctr)) + 1)),
                          DEPva = apply(expmtr, 1, protein_de_via_wilcox))
      finalres <- rbind(finalres, deres)
    }
    
  }
}

save(finalres, file = "../8.results/tf.regu.sub.RData")

load("../8.results/tf.regu.sub.RData")

finalres <- finalres[finalres$Meantest > 1, ]
finalres <- finalres[finalres$DEFC > 0, ]
finalres <- finalres[finalres$DEPva < 0.01, ]

tsggenes <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F, sep = "\t")
finalres <- finalres[finalres$ID %in% tsggenes$Gene, ]



tsggenes <- fread("../0.dataset/Human_TSGs.txt", header = T, stringsAsFactors = F, data.table = F, sep = "\t")
nucmut <- fread("../4.panMut/singmut.delta.sig.annote.Nucleus", 
                header = T, stringsAsFactors = F, data.table = F)
nucmut <- nucmut[nucmut$gene %in% tsggenes$GeneSymbol, ]
nucmut <- nucmut[nucmut$deltaScore < -0.1, ]
alltfs <- fread("../0.dataset/TF-Target-information.txt", header = T, stringsAsFactors = F, data.table = F)
nucmut <- nucmut[nucmut$gene %in% alltfs$TF, ]
tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)

tfs <- "E2F1"
tumors = "STAD"
allmuts <- fread(paste("/home/yukai/data/TCGA/mutect_tcga/TCGA-", tumors, ".mutect2_snv.tsv", sep = ""), 
                 header = T, stringsAsFactors = F, data.table = F)
mutpairs <- nucmut[nucmut$gene == tfs, ]$Amino_Acid_Change
subgenes <- unique(alltfs[alltfs$TF == tfs, ]$target)
mutsams <- unique(allmuts[allmuts$Amino_Acid_Change %in% mutpairs, ]$Sample_ID)
nonsams <- unique(allmuts[!(allmuts$Amino_Acid_Change %in% mutpairs), ]$Sample_ID)
tmp_clininfo <- data.frame(ID = c(mutsams, nonsams),
                           types = c(rep("Mut", length(mutsams)), rep("Non", length(nonsams))))

expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", tumors, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]
expmtr <- expmtr[rownames(expmtr) %in% subgenes, ]
exp_tes <- expmtr[, names(expmtr) %in% mutsams]
exp_ctr <- expmtr[, names(expmtr) %in% nonsams]
res4plot1 <- data.frame(ID = c(names(expmtr)),
                       Pros = c(rep("LYN", ncol(expmtr))),
                       Class = c(rep("E2F1", ncol(expmtr))),
                       ExpValue = c(as.numeric(expmtr["LYN", ])))
res4plot1 <- merge(res4plot1, tmp_clininfo, by.x = "ID", by.y = "ID")

tfs <- "E2F1"
tumors = "COAD"
allmuts <- fread(paste("/home/yukai/data/TCGA/mutect_tcga/TCGA-", tumors, ".mutect2_snv.tsv", sep = ""), 
                 header = T, stringsAsFactors = F, data.table = F)
mutpairs <- nucmut[nucmut$gene == tfs, ]$Amino_Acid_Change
subgenes <- unique(alltfs[alltfs$TF == tfs, ]$target)
mutsams <- unique(allmuts[allmuts$Amino_Acid_Change %in% mutpairs, ]$Sample_ID)
nonsams <- unique(allmuts[!(allmuts$Amino_Acid_Change %in% mutpairs), ]$Sample_ID)
tmp_clininfo <- data.frame(ID = c(mutsams, nonsams),
                           types = c(rep("Mut", length(mutsams)), rep("Non", length(nonsams))))

expmtr <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-", tumors, ".htseq_fpkm.tsv.cv.txt", sep = ""), 
                header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$Ensembl_ID
expmtr <- expmtr[, -1]
expmtr <- expmtr[rownames(expmtr) %in% subgenes, ]
exp_tes <- expmtr[, names(expmtr) %in% mutsams]
exp_ctr <- expmtr[, names(expmtr) %in% nonsams]
res4plot2 <- data.frame(ID = c(names(expmtr)),
                       Pros = c(rep("TRIM7", ncol(expmtr))),
                       Class = c(rep("E2F1", ncol(expmtr))),
                       ExpValue = c(as.numeric(expmtr["TRIM7", ])))
res4plot2 <- merge(res4plot2, tmp_clininfo, by.x = "ID", by.y = "ID")

res4plot <- rbind(res4plot1, res4plot2)

p1 <- ggplot(res4plot, aes(x = types, y = ExpValue, fill = types))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", method = "wilcox")+
  ylab("mRNA expression level")+
  theme_classic2()
q <- facet(p1, facet.by = "Pros", scales = "free_y")
q
pdf("../8.results/01.tfs.dn.target.pdf", width = 7, height =4)
print(q)
dev.off()

#################################Mutation effect selection#################################
########affcted mutations########
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
signuc <- genesmerge[genesmerge$Class == "Nucleus", ]
sigcyt <- genesmerge[genesmerge$Class == "Cytosol", ]
sigmem <- genesmerge[genesmerge$Class == "Membrane", ]
sigmit <- genesmerge[genesmerge$Class == "Mitochondrion" & abs(genesmerge$Delta) > 0.8, ]
sigsec <- genesmerge[genesmerge$Class == "Secreted" & abs(genesmerge$Delta) > 0.6, ]
hsdataset <- fread("../4.panMut/hs.protein.score.Cytosol", header = F, stringsAsFactors = F, data.table = F)

res <- data.frame()

sigregions <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
singmut <- merge(signuc, hsdataset, by.x = "gene", by.y = "V2")
commonpros <- intersect(sigregions$V1, singmut$V1)
genes <- c()
mutreg <- c()
mutoth <- c()
num = 1
for (pro in commonpros) {
  tmpmut <- singmut[singmut$V1 == pro, ]
  mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(gsub(".*_", "", tmpmut$ID), "[^0-9]+")))))
  proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  otherregion <- setdiff(proseq, nlsregion)
  nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
  otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
  genes <- c(genes, pro)
  mutreg <- c(mutreg, nlsrate)
  mutoth <- c(mutoth, otherrate)
  print(paste(num, " of ", length(commonpros), " has been finished!!!"))
  num = num + 1
}
tmpres <- data.frame(ID = genes,
                     MutNLS = mutreg,
                     MutOth = mutoth,
                     Class = "Nucleus",
                     stringsAsFactors = F)
res <- rbind(res, tmpres)


sigregions <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
singmut <- merge(sigmit, hsdataset, by.x = "gene", by.y = "V2")
commonpros <- intersect(sigregions$V1, singmut$V1)
genes <- c()
mutreg <- c()
mutoth <- c()
num = 1
for (pro in commonpros) {
  tmpmut <- singmut[singmut$V1 == pro, ]
  mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(gsub(".*_", "", tmpmut$ID), "[^0-9]+")))))
  proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  otherregion <- setdiff(proseq, nlsregion)
  nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
  otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
  genes <- c(genes, pro)
  mutreg <- c(mutreg, nlsrate)
  mutoth <- c(mutoth, otherrate)
  print(paste(num, " of ", length(commonpros), " has been finished!!!"))
  num = num + 1
}
tmpres <- data.frame(ID = genes,
                     MutNLS = mutreg,
                     MutOth = mutoth,
                     Class = "Mito",
                     stringsAsFactors = F)
res <- rbind(res, tmpres)

sigregions <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
singmut <- merge(sigsec, hsdataset, by.x = "gene", by.y = "V2")
commonpros <- intersect(sigregions$V1, singmut$V1)
genes <- c()
mutreg <- c()
mutoth <- c()
num = 1
for (pro in commonpros) {
  tmpmut <- singmut[singmut$V1 == pro, ]
  mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(gsub(".*_", "", tmpmut$ID), "[^0-9]+")))))
  proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
  otherregion <- setdiff(proseq, nlsregion)
  nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
  otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
  genes <- c(genes, pro)
  mutreg <- c(mutreg, nlsrate)
  mutoth <- c(mutoth, otherrate)
  print(paste(num, " of ", length(commonpros), " has been finished!!!"))
  num = num + 1
}
tmpres <- data.frame(ID = genes,
                     MutNLS = mutreg,
                     MutOth = mutoth,
                     Class = "Secreted",
                     stringsAsFactors = F)
res <- rbind(res, tmpres)

res4plot <- melt(res, id.vars = c("ID", "Class"))
res4plot$Class <- factor(res4plot$Class, c("Nucleus", "Mito", "Secreted"))
p4 <- ggplot(res4plot, aes(x = Class, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation counts per 1000 AAs")+
  theme_linedraw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p4

ggsave("../8.results/02.mut.delta.region.pdf", p4, width = 6, height = 3.5)


######## mut freq and survival########
allppi <- fread("/home/yukai/data/exp_ppi/exp_uniq.txt", header = F)
####mut freq####
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
signuc <- genesmerge[genesmerge$Class == "Nucleus", ]

tumorlist <- fread("../0.dataset/select_cancer.list", header = F, stringsAsFactors = F, data.table = F)
res <- data.frame()
for (cancer in tumorlist$V1) {
  singsurv <- fread(paste("/home/yukai/data/Datasets/TCGA20201022/survival/TCGA-", cancer, ".survival.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singsurv <- singsurv[as.numeric(substr(singsurv$sample, 14, 15)) < 10, ]
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  sigregions <- signuc[signuc$Delta > 0, ]
  commonpros <- intersect(sigregions$gene, singmut$gene)
  genes <- c()
  mutsite <- c()
  mutsam <- c()
  num = 1
  for (pro in commonpros) {
    tmpregion <- sigregions[sigregions$gene == pro, ]
    tmpmut <- singmut
    tmpmut$ID <- paste(tmpmut$gene, tmpmut$Amino_Acid_Change, sep = "_")
    newtmpmut <- tmpmut[tmpmut$Amino_Acid_Change %in% gsub(".*_", "", sigregions$ID), ]
    genes <- c(genes, pro)
    mutsite <- c(mutsite, length(unique(newtmpmut$Sample_ID)))
    mutsam <- c(mutsam, nrow(singsurv))
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutSite = mutsite,
                       MutSam = mutsam,
                       Tumor = cancer,
                       Class = "Increased",
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
  sigregions <- signuc[signuc$Delta < 0, ]
  commonpros <- intersect(sigregions$gene, singmut$gene)
  genes <- c()
  mutsite <- c()
  mutsam <- c()
  num = 1
  for (pro in commonpros) {
    tmpregion <- sigregions[sigregions$gene == pro, ]
    tmpmut <- singmut
    tmpmut$ID <- paste(tmpmut$gene, tmpmut$Amino_Acid_Change, sep = "_")
    newtmpmut <- tmpmut[tmpmut$Amino_Acid_Change %in% gsub(".*_", "", tmpregion$ID), ]
    genes <- c(genes, pro)
    mutsite <- c(mutsite, length(unique(newtmpmut$Sample_ID)))
    mutsam <- c(mutsam, nrow(singsurv))
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutSite = mutsite,
                       MutSam = mutsam,
                       Tumor = cancer,
                       Class = "Decreased",
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}

write.table(res, "../8.results/02.mut.nuc.gene.count.txt", row.names = F, col.names = T, quote = F, sep = "\t")

####mut survival####
mutsurv <- fread("../8.results/02.mut.nuc.gene.count.txt", header = T)
mutsurv <- mutsurv[mutsurv$MutSite < 31 & mutsurv$MutSite > 9, ]
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
signuc <- genesmerge[genesmerge$Class == "Nucleus", ]

hrs <- c()
pvas <- c()
for (i in c(1:nrow(mutsurv))) {
  selgene <- mutsurv[i, ]$ID
  selcancer <- mutsurv[i, ]$Tumor
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", selcancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  allsamples <- unique(singmut$Sample_ID)
  if (mutsurv[i, ]$Class == "Increased") {
    sigregions <- signuc[signuc$Delta > 0, ]
  }else{
    sigregions <- signuc[signuc$Delta < 0, ]
  }
  tmpregion <- sigregions[sigregions$gene == selgene, ]
  tmpmut <- singmut
  tmpmut$ID <- paste(tmpmut$gene, tmpmut$Amino_Acid_Change, sep = "_")
  newtmpmut <- tmpmut[tmpmut$Amino_Acid_Change %in% gsub(".*_", "", tmpregion$ID), ]
  mutsamples <- unique(newtmpmut$Sample_ID)
  survinfo <- fread(paste("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/survival/TCGA-", selcancer, ".survival.tsv", sep = ""), 
                    header = T, stringsAsFactors = F, data.table = F)
  survinfo <- survinfo[as.numeric(substr(survinfo$sample, 14, 15)) < 10, ]
  rownames(survinfo) <- survinfo$sample
  commonsams <- intersect(rownames(survinfo), allsamples)
  survinfo <- survinfo[commonsams, ]
  cluster_surv <- survinfo
  cluster_surv$Type = ifelse(cluster_surv$sample %in% mutsamples, "1.Mut", "0.Others")
  tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
  hrs <- c(hrs, tmp$conf.int[[1]])
  pvas <- c(pvas, tmp$logtest[[3]])
}
mutsurv$HR <- hrs
mutsurv$Pval <- pvas

write.table(mutsurv, "../8.results/02.mut.nuc.gene.surv.txt", row.names = F, col.names = T, quote = F, sep = "\t")


tmpsurv <- mutsurv[mutsurv$HR > 1 & mutsurv$Pval < 0.05, ]

########nuclear select########
allppi <- fread("/home/yukai/data/exp_ppi/exp_uniq.txt", header = F)
alltsgs <- fread("../0.dataset/Human_TSGs.txt", header = T, stringsAsFactors = F)

genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)
signuc <- genesmerge[genesmerge$Class == "Nucleus", ]

nuc_n2y <- signuc[signuc$Delta > 0.1, ]
nuc_y2n <- signuc[signuc$Delta < -0.1, ]

##TSG non-nuc to nuc
genes <- intersect(alltsgs$GeneSymbol, setdiff(nuc_n2y$gene, nuc_y2n$gene))

alldata <- fread("../../ProNuclear/3.model/roc_data.alldata.txt", header =T, stringsAsFactors = F, data.table = F)
cutoff=alldata$score[which.max(alldata$sn-1+alldata$sp)]
res <- data.frame(ID = genes)

nucloc <- fread("../4.panMut/hs.protein.score.Nucleus", header = F)
nucloc$Loc <- ifelse(nucloc$V4 > cutoff, "Nucleus", "Others")
res <- merge(res, nucloc[, c(2, 5)], by.x = "ID", by.y = "V2")

res <- res[res$Loc == "Others", ]

resppi1 <- merge(res, allppi[, c(3,4)], by.x = "ID", by.y = "V3")
resppi2 <- merge(res, allppi[, c(3,4)], by.x = "ID", by.y = "V4")
names(resppi1) <- c("ID", "IDLoc", "PPIs")
names(resppi2) <- c("ID", "IDLoc", "PPIs")
resppi <- rbind(resppi1, resppi2)
resppi <- unique(resppi)

resppiloc <- merge(resppi, nucloc[, c(2, 5)], by.x = "PPIs", by.y = "V2")

##enrichment
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_description, gene_symbol)
ewp1 <- enricher(resppiloc[resppiloc$Loc == "Nucleus", ]$PPIs, pvalueCutoff = 1,
                 TERM2GENE = kegg, qvalueCutoff = 1, minGSSize = 1)
ewp2 <- enricher(resppiloc[resppiloc$Loc == "Others", ]$PPIs, pvalueCutoff = 1,
                 TERM2GENE = kegg, qvalueCutoff = 1, minGSSize = 1)
sels1 <- c(1,3,7,9)
sels2 <- c(2,6,8,9)
allpath <- data.frame(ID = c(ewp1@result[sels1, ]$ID, ewp2@result[sels2, ]$ID),
                      GR = c(as.numeric(gsub("/.*", "", ewp1@result[sels1, ]$GeneRatio)) / as.numeric(gsub(".*/", "", ewp1@result[sels1, ]$GeneRatio)), 
                             - as.numeric(gsub("/.*", "", ewp2@result[sels2, ]$GeneRatio)) / as.numeric(gsub(".*/", "", ewp2@result[sels2, ]$GeneRatio))),
                      pva = c(-log10(ewp1@result[sels1, ]$p.adjust), -log10(ewp2@result[sels2, ]$p.adjust)),
                      stringsAsFactors = F)

allpath <- allpath[order(allpath$GR), ]
allpath$ID <- factor(allpath$ID, allpath$ID)
p1 <- ggplot(allpath, aes(x = GR, y = ID, fill = pva))+
  geom_col(width = 0.8)+
  #stat_compare_means(label = "p.signif")+
  ylab("Relative Rates")+
  theme_classic2()
p1
ggsave("../8.results/02.mut.ppi.enrich.pdf", p1, width = 6, height = 3)

##annotation (hub genes, CAGs, TSGs, nuclear, each select 4 pathways)
ppis4net <- resppiloc
ppis4net <- ppis4net[ppis4net$ID != ppis4net$PPIs, ]
#hub genes
anno1 <- data.frame(ID = unique(ppis4net$ID),
                    Class = "Hub Genes")
#other genes
allppis <- data.frame(ID = unique(ppis4net$PPIs))

cags <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F)
tsgs <- fread("../0.dataset/Human_TSGs.txt", header = T, stringsAsFactors = F)
nucle <- ppis4net
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_description, gene_symbol)
kegg <- kegg[kegg$gs_description %in% allpath$ID, ]
names(kegg) <- c("Pathways", "Genes")

allppis$CAGs <- ifelse(allppis$ID %in% cags$Gene, "CAGs", "")
allppis$TSGs <- ifelse(allppis$ID %in% tsgs$GeneSymbol, "TSGs", "")
allppis$Nucs <- ifelse(allppis$ID %in% nucle[nucle$Loc == "Nucleus", ]$PPIs, "Nucleus", "")
allppis <- merge(allppis, kegg, by.x = "ID", by.y = "Genes", all.x = T)

anno2 <- melt(allppis, id.vars = "ID")
anno2 <- na.omit(anno2[, c(1, 3)])
anno2 <- anno2[anno2$value != "", ]
names(anno2) <- c("ID", "Class")
anno2 <- unique(anno2)
anno <- rbind(anno1, anno2)

write.table(anno, "../8.results/02.mut.ppi.net.annot.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(resppiloc, "../8.results/02.mut.ppi.net.ppi.txt", row.names = F, col.names = T, quote = F, sep = "\t")


##TNK1
resppiloc <- fread("../8.results/02.mut.ppi.net.ppi.txt", header = T, stringsAsFactors = F, data.table = F)
allppis <- resppiloc
allppis <- allppis[allppis$ID != allppis$PPIs, ]
selppis <- allppis[allppis$ID == "TNK1", ]

allpros <- fread("../8.results/02.mut.ppi.net.annot.txt", header = T)
selpros <- unique(c(selppis$ID, selppis$PPIs))
othpros <- setdiff(allpros$ID, selpros)

write.table(selpros, "../8.results/02.mut.ppi.selnet.tnk1.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(othpros, "../8.results/02.mut.ppi.selnet.nontnk1.txt", row.names = F, col.names = F, quote = F, sep = "\t")


ppis4net <- resppiloc
ppis4net <- ppis4net[ppis4net$ID != ppis4net$PPIs, ]
allppis <- data.frame(ID = unique(ppis4net$PPIs))
cags <- fread("/home/yukai/data/Datasets/HPA/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F)
tsgs <- fread("../0.dataset/Human_TSGs.txt", header = T, stringsAsFactors = F)
nucle <- ppis4net
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_description, gene_symbol)
kegg <- kegg[kegg$gs_description %in% allpath$ID, ]
names(kegg) <- c("Pathways", "Genes")
allppis$CAGs <- ifelse(allppis$ID %in% cags$Gene, "CAGs", "")
allppis$TSGs <- ifelse(allppis$ID %in% tsgs$GeneSymbol, "TSGs", "")
allppis$Nucs <- ifelse(allppis$ID %in% nucle[nucle$Loc == "Nucleus", ]$PPIs, "Nucleus", "")
allppis <- merge(allppis, kegg, by.x = "ID", by.y = "Genes", all.x = T)

selppis <- merge(selppis, allppis, by.x = "PPIs", by.y = "ID")







#################################RNAseq results#################################
########PCA########
expmtr <- fread("../9.seqData/2.MergeQuant/Count_symbol.txt", header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$V1
expmtr <- expmtr[, -1]
tmpmtr <- expmtr
df_pca <- prcomp(tmpmtr)
pcaData_cgs <- as.data.frame(df_pca$rotation)
pcaData_cgs$label <- rownames(pcaData_cgs)
df_out <- as.data.frame(df_pca$x)
percentage_cgs <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage_cgs <- paste( colnames(df_out), "(", paste( as.character(percentage_cgs), "%", ")", sep="") )
p3 <- ggplot(pcaData_cgs, aes(PC1, PC2, color=gsub("-.", "", rownames(pcaData_cgs)))) +
  geom_point(size=3) +
  stat_ellipse(level = 0.95, show.legend = F)+
  ggtitle("RNAseq")+
  #geom_mark_hull()+
  xlab(percentage_cgs[1]) + ylab(percentage_cgs[2])+
  scale_colour_hue("Type") +
  theme_bw()+
  theme(plot.title=element_text(size = '16', color = 'black', hjust = 0.5))
p3
p4 <- ggplot(pcaData_cgs, aes(PC1, PC2, color=gsub("-.", "", rownames(pcaData_cgs)))) +
  geom_point(size=3) +
  geom_label(aes(label = label))+
  stat_ellipse(level = 0.95, show.legend = F)+
  ggtitle("RNAseq")+
  #geom_mark_hull()+
  xlab(percentage_cgs[1]) + ylab(percentage_cgs[2])+
  scale_colour_hue("Type") +
  theme_bw()+
  theme(plot.title=element_text(size = '16', color = 'black', hjust = 0.5))
p4

skipsamples <- c("DLD1R14M-1", "DLD1WT-3", "HCT116R14M-1", "HCT116WT-1")
tmpmtr <- expmtr[, !(names(expmtr) %in% skipsamples)]
df_pca <- prcomp(tmpmtr)
pcaData_cgs <- as.data.frame(df_pca$rotation)
pcaData_cgs$label <- rownames(pcaData_cgs)
df_out <- as.data.frame(df_pca$x)
percentage_cgs <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage_cgs <- paste( colnames(df_out), "(", paste( as.character(percentage_cgs), "%", ")", sep="") )
p5 <- ggplot(pcaData_cgs, aes(PC1, PC2, color=gsub("-.", "", rownames(pcaData_cgs)))) +
  geom_point(size=3) +
  stat_ellipse(level = 0.95, show.legend = F)+
  ggtitle("RNAseq")+
  #geom_mark_hull()+
  xlab(percentage_cgs[1]) + ylab(percentage_cgs[2])+
  scale_colour_hue("Type") +
  theme_bw()+
  theme(plot.title=element_text(size = '16', color = 'black', hjust = 0.5))
p5
ggsave("../9.seqData/3.Deseq/rna.pca.pdf", p3, width = 5, height = 4)
ggsave("../9.seqData/3.Deseq/rna.pca.label.pdf", p4, width = 5, height = 4)
ggsave("../9.seqData/3.Deseq/rna.pca.skip.pdf", p5, width = 5, height = 4)
########deanalysis########
expmtr <- fread("../9.seqData/2.MergeQuant/Count_symbol.txt", header = T, stringsAsFactors = F, data.table = F)
rownames(expmtr) <- expmtr$V1
expmtr <- expmtr[, -1]

skipsamples <- c("DLD1R14M-1", "DLD1WT-3", "HCT116R14M-1", "HCT116WT-1")
##DLD1
desmtr <- expmtr[, 1:6]
desmtr <- desmtr[, !(names(desmtr) %in% skipsamples)]
library(DESeq2)
colData <- data.frame(types = gsub("DLD1", "", gsub("-.", "", names(desmtr))))
rownames(colData) <- names(desmtr)
countData = round(desmtr)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ types)
dds <- DESeq(dds)
type_level <- levels(as.factor(colData$types))
comb <- combn(type_level,2)
res <- results(dds, contrast=c("types","R14M","WT"))
res <- as.data.frame(res)
res <- na.omit(res)
write.table(res, file = "../9.seqData/3.Deseq/rna.dld1.deseq.xls",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)
##HCT116
desmtr <- expmtr[, 7:12]
desmtr <- desmtr[, !(names(desmtr) %in% skipsamples)]
library(DESeq2)
colData <- data.frame(types = gsub("HCT116", "", gsub("-.", "", names(desmtr))))
rownames(colData) <- names(desmtr)
countData = round(desmtr)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ types)
dds <- DESeq(dds)
type_level <- levels(as.factor(colData$types))
comb <- combn(type_level,2)
res <- results(dds, contrast=c("types","R14M","WT"))
res <- as.data.frame(res)
res <- na.omit(res)
write.table(res, file = "../9.seqData/3.Deseq/rna.hct116.deseq.xls",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

########venn########
genes <- fread("/home/yukai/data/hg38/EnsID2Symbol.txt", header = F)
allgenes <- genes[genes$V3 == "protein_coding", ]$V2

de_dld1 <- fread("../9.seqData/3.Deseq/rna.dld1.deseq.xls", header = T, stringsAsFactors = F, data.table = F)
de_hct116 <- fread("../9.seqData/3.Deseq/rna.hct116.deseq.xls", header = T, stringsAsFactors = F, data.table = F)

de_dld1 <- de_dld1[de_dld1$V1 %in% allgenes, ]
de_hct116 <- de_hct116[de_hct116$V1 %in% allgenes, ]

lfc <- log(2)
pva <- 0.05

up_dld1 <- de_dld1[de_dld1$log2FoldChange > lfc & de_dld1$padj < pva, ]$V1
dn_dld1 <- de_dld1[de_dld1$log2FoldChange < -lfc & de_dld1$padj < pva, ]$V1
up_hct116 <- de_hct116[de_hct116$log2FoldChange > lfc & de_hct116$padj < pva, ]$V1
dn_hct116 <- de_hct116[de_hct116$log2FoldChange < -lfc & de_hct116$padj < pva, ]$V1

library(VennDiagram)
venn.plot <- venn.diagram(list(DLD1_Up=up_dld1,
                               DLD1_Dn=dn_dld1,
                               HCT_Up=up_hct116,
                               HCT_Dn=dn_hct116),
                          resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5,0.5),
                          #fill=c("red", "blue"),
                          main="",
                          main.cex = 2, 
                          filename = NULL)

pdf(file="../9.seqData/3.Deseq/dld1.hct.de.overlap.pdf")
grid.draw(venn.plot)
dev.off()




########pathways########
genes <- fread("/home/yukai/data/hg38/EnsID2Symbol.txt", header = F)
allgenes <- genes[genes$V3 == "protein_coding", ]$V2

de_dld1 <- fread("../9.seqData/3.Deseq/rna.dld1.deseq.xls", header = T, stringsAsFactors = F, data.table = F)
de_dld1 <- de_dld1[de_dld1$V1 %in% allgenes, ]
lfc <- log(2)
pva <- 0.05
up_dld1 <- de_dld1[de_dld1$log2FoldChange > lfc & de_dld1$padj < pva, ]$V1
dn_dld1 <- de_dld1[de_dld1$log2FoldChange < -lfc & de_dld1$padj < pva, ]$V1

library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_description, gene_symbol)

##correlation GSEA
kegg_up <- enricher(up_dld1, TERM2GENE = kegg, pvalueCutoff = 1,minGSSize = 1)
write.table(kegg_up@result, file = '../9.seqData/4.Pathway/path.kegg.up.txt',
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)
kegg_dn <- enricher(dn_dld1, TERM2GENE = kegg, pvalueCutoff = 1,minGSSize = 1)
write.table(kegg_dn@result, file = '../9.seqData/4.Pathway/path.kegg.dn.txt',
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

##sep
ewp1 <- fread("../9.seqData/4.Pathway/path.kegg.up.txt", header = T, stringsAsFactors = F, data.table = F)
ewp1$ID <- paste(ewp1$ID, ".up", sep = "")
ewp2 <- fread("../9.seqData/4.Pathway/path.kegg.dn.txt", header = T, stringsAsFactors = F, data.table = F)
ewp2$ID <- paste(ewp2$ID, ".dn", sep = "")
sels1 <- c(4,6,8,10,11)
sels2 <- c()
allpath <- data.frame(ID = c(ewp1[sels1, ]$ID, ewp2[sels2, ]$ID),
                      GR = c(rep("0.up", length(sels1)), rep("1.dn", length(sels2))),
                      pva = c(-log10(ewp1[sels1, ]$p.adjust), log10(ewp2[sels2, ]$p.adjust)),
                      stringsAsFactors = F)

allpath <- allpath[order(allpath$pva), ]
allpath$ID <- factor(allpath$ID, allpath$ID)
p1 <- ggplot(allpath, aes(x = pva, y = ID, fill = GR))+
  geom_col(width = 0.8)+
  #stat_compare_means(label = "p.signif")+
  ylab("Relative Rates")+
  theme_classic2()
p1
ggsave("../9.seqData/4.Pathway/path.overlap.enrich.pdf", p1, width = 6, height = 3)

########heatmap########
genes <- fread("/home/yukai/data/hg38/EnsID2Symbol.txt", header = F)
allgenes <- genes[genes$V3 == "protein_coding", ]$V2

ewp1 <- fread("../9.seqData/4.Pathway/path.kegg.up.txt", header = T, stringsAsFactors = F, data.table = F)
selgenes <- strsplit(ewp1[10, ]$geneID, "/")[[1]]
res <- fread("../9.seqData/3.Deseq/rna.dld1.deseq.xls", header = T, stringsAsFactors = F, data.table = F)
res <- res[res$V1 %in% allgenes, ]
desmtr <- fread("../9.seqData/2.MergeQuant/TPM_symbol.txt", header = T, data.table = F)
skipsamples <- c("DLD1R14M-1", "DLD1WT-3", "HCT116R14M-1", "HCT116WT-1")
##DLD1
desmtr <- expmtr[, 1:6]
desmtr <- desmtr[, !(names(desmtr) %in% skipsamples)]

fc <- 1.5
lfc <- log2(fc)
pval <- 0.05

sigGenes_up = res[(res$log2FoldChange > lfc & res$padj < pval), ]$V1
sigGenes_down = res[(res$log2FoldChange < -lfc & res$padj < pval), ]$V1

library(ComplexHeatmap)
demtr <- res[res$V1 %in% c(sigGenes_up, sigGenes_down), ]
demtr <- demtr[order(demtr$log2FoldChange, decreasing = T), ]
mtr4heat <- desmtr[demtr$V1, ]
mark <- selgenes
colData <- data.frame(row.names = colnames(mtr4heat),
                      types = gsub("DLD1", "", gsub("-.", "", names(desmtr))))
group_col <- c("#bfb2d5","#f1937f")
names(group_col) <- c('R14M','WT')
col <- list(types = group_col)

ht = HeatmapAnnotation(
  Subtype = colData$types,
  col = list(Subtype = c("WT" = "#bfb2d5", "R14M" = "#f1937f")
  ))
ha = rowAnnotation(foo = anno_mark(at = match(mark,rownames(mtr4heat)), 
                                   labels = mark))
scaled_mat = t(scale(t(mtr4heat)))
Heatmap(scaled_mat, name = "mat", cluster_rows = FALSE, 
        top_annotation = ht, right_annotation = ha,
        col = colorRampPalette(c("#06a7cd", "white", "#e74a32"))(100),
        show_row_names = F, show_column_names = F)
pdf("../9.seqData/4.Pathway/degenes.heatmap.pdf", width = 5, height = 4)
Heatmap(scaled_mat, name = "mat", cluster_rows = FALSE, 
        top_annotation = ht, right_annotation = ha,
        col = colorRampPalette(c("#06a7cd", "white", "#e74a32"))(100),
        show_row_names = F, show_column_names = F)
dev.off()




#################################10.ReNuclear#################################
########mutation enriched in signals for combined########
sigregions1 <- fread("../0.dataset/valid.nls.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions2 <- fread("../0.dataset/valid.nes.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions3 <- fread("../0.dataset/sp.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions4 <- fread("../0.dataset/mito.human.pos.txt", header = F, stringsAsFactors = F, data.table = F)
sigregions <- rbind(sigregions1, sigregions2)
pro2seq <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset <- fread("../4.panMut/hs.protein.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
hsdataset$GeneName <- hsdataset$V2
tumorlist <- c("Pancancer")
res <- data.frame()
for (cancer in tumorlist) {
  singmut <- fread(paste("../../ProNuclear/6.pancancer/1.mutfile/TCGA-", cancer, ".mutect2_snv.tsv", sep = ""), 
                   header = T, stringsAsFactors = F, data.table = F)
  singmut <- merge(singmut, hsdataset[, c(1,5)], by.x = "gene", by.y = "GeneName")
  commonpros <- intersect(sigregions$V1, singmut$V1)
  genes <- c()
  mutreg <- c()
  mutoth <- c()
  num = 1
  for (pro in commonpros) {
    tmpmut <- singmut[singmut$V1 == pro, ]
    mutpos <- unique(na.omit(as.numeric(unlist(strsplit2(tmpmut$Amino_Acid_Change, "[^0-9]+")))))
    proseq <- c(1:nchar(hsdataset[hsdataset$V1 == pro, ]$V3))
    nlsregion <- str2lst(sigregions[sigregions$V1 == pro, ]$V2, nchar(pro2seq[pro2seq$V1 == pro, ]$V3))
    otherregion <- setdiff(proseq, nlsregion)
    nlsrate <- log(((length(intersect(mutpos, nlsregion))+1)/length(nlsregion))*1000)
    otherrate <- log(((length(intersect(mutpos, otherregion))+1)/length(otherregion))*1000)
    genes <- c(genes, pro)
    mutreg <- c(mutreg, nlsrate)
    mutoth <- c(mutoth, otherrate)
    print(paste(num, " of ", length(commonpros), " has been finished!!!"))
    num = num + 1
  }
  
  tmpres <- data.frame(ID = genes,
                       MutNLS = mutreg,
                       MutOth = mutoth,
                       Tumor = cancer,
                       stringsAsFactors = F)
  res <- rbind(res, tmpres)
}
write.table(res, "../10.NucleRes/inuloc.point.mut.Mergedregion.txt", col.names = T, row.names = F, sep = '\t', quote = F)

res <- fread("../10.NucleRes/inuloc.point.mut.Mergedregion.txt", header = T, stringsAsFactors = F, data.table = F)
res$MutNLS <- log10(exp(res$MutNLS))
res$MutOth <- log10(exp(res$MutOth))
res4plot4 <- melt(res, id.vars = c("ID", "Tumor"))
p1 <- ggplot(res4plot4, aes(x = variable, y = value, fill = variable))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif")+
  ylab("log10 Mutation count per 1000 AAs (MT)")+
  theme_classic2()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("../10.NucleRes/inuloc.point.mut.distibute.Merged.boxplot.pdf", width = 4, height = 4)
print(p1)
dev.off()

########sample distribution of locations########
selorgas <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Saccharomyces cerevisiae", "Drosophila melanogaster")
usedlocs <- c("Nucleus")
#colorinfo <- brewer.pal(5, "Set1")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")
subloc <- subloc[subloc$V3 %in% usedlocs, ]
subloc <- subloc[subloc$V2 %in% selorgas, ]

subloc$V2 <- factor(subloc$V2, selorgas)
p <- ggplot(subloc, aes(V2, fill=V3))+
  geom_bar(position = "stack")+
  #scale_fill_manual(values = colorinfo)+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p
pdf("../10.NucleRes/inuloc.sublocs.distribut.multiloc.organism.pdf", width = 5, height =4)
print(p)
dev.off()


########amino acid composition for all locs########
selorgas <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Saccharomyces cerevisiae", "Drosophila melanogaster")
usedlocs <- c("Nucleus", "Cytosol", "Membrane", "Mitochondrion", "Secreted")
AA_list_k=c("H","R","K","F","A","L","M","I","W","P","V","C","G","Q","N","S","Y","T","D","E")
colorinfo <- brewer.pal(5, "Set1")
subloc <- fread("../0.data/0.1uniprot.subloc.model_organism.txt", header = F, stringsAsFactors = F, data.table = F)
subloc <- subloc %>% separate_rows(V3, sep = ",")

suball <- subloc
subnuc <- subloc[subloc$V3 == "Nucleus", ]
subcyt <- subloc[subloc$V3 == "Cytosol", ]
submem <- subloc[subloc$V3 == "Membrane", ]
submit <- subloc[subloc$V3 == "Mitochondrion", ]
subsec <- subloc[subloc$V3 == "Secreted", ]

suballres <- as.data.frame(table(strsplit(paste(suball$V4, collapse = ""), "")[[1]]))
suballres$Rate <- suballres$Freq / sum(suballres$Freq)
subnucres <- as.data.frame(table(strsplit(paste(subnuc$V4, collapse = ""), "")[[1]]))
subnucres$Rate <- subnucres$Freq / sum(subnucres$Freq)
subcytres <- as.data.frame(table(strsplit(paste(subcyt$V4, collapse = ""), "")[[1]]))
subcytres$Rate <- subcytres$Freq / sum(subcytres$Freq)
submemres <- as.data.frame(table(strsplit(paste(submem$V4, collapse = ""), "")[[1]]))
submemres$Rate <- submemres$Freq / sum(submemres$Freq)
submitres <- as.data.frame(table(strsplit(paste(submit$V4, collapse = ""), "")[[1]]))
submitres$Rate <- submitres$Freq / sum(submitres$Freq)
subsecres <- as.data.frame(table(strsplit(paste(subsec$V4, collapse = ""), "")[[1]]))
subsecres$Rate <- subsecres$Freq / sum(subsecres$Freq)

res <- rbind(suballres, subnucres)
res$Class <- c(rep("BG", nrow(suballres)), rep("Nuc", nrow(subnucres)))
res <- res[res$Var1 %in% AA_list_k, ]
res$Var1 <- factor(res$Var1, AA_list_k)
p6 <- ggplot(res, aes(x = Var1, y = Rate, group = Class))+
  geom_line(color = "gray50", alpha = 0.5)+
  geom_point(aes(shape=Class,colour=Var1), size = 5)+
  theme_classic2()
p6

pdf("../10.NucleRes/inuloc.sublocs.distribut.aa.freq.pdf", width = 10, height =7)
print(p6)
dev.off()




########SYSUCC_CRC_Mutation########
allmuts <- fread("../10.NucleRes/SYSUCC_CRC_Mutation_withSilent.maf",
                 header = T, stringsAsFactors = F, data.table = F)
allmuts <- allmuts[allmuts$Variant_Classification == "Missense_Mutation", ]
allmuts <- allmuts[allmuts$Variant_Classification == "Missense_Mutation", ]
newres <- unique(allmuts[, c("Hugo_Symbol", "HGVSp_Short", "Tumor_Sample_Barcode")])
newres$ID <- paste(newres$Hugo_Symbol, newres$HGVSp_Short, sep = "_")

allsams <- length(unique(allmuts$Tumor_Sample_Barcode))
rescount <- as.data.frame(table(newres$ID))
rescount$Ratio <- (rescount$Freq / allsams) * 100

rescount$Hugo_Symbol <- strsplit2(rescount$Var1, "_")[, 1]
rescount$HGVSp_Short <- strsplit2(rescount$Var1, "_")[, 2]

selmuts <- c("KAT8", "DEAF1", "NFE2L1", "CMAS", "NHLRC1", "PTEN", "CHFR")
selrescount <- rescount[rescount$Hugo_Symbol %in% selmuts, ]
write.table(selrescount, "../10.NucleRes/sysucc.selrescount.uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)


rescount <- rescount[rescount$Freq > 1, ]
write.table(rescount, "../10.NucleRes/sysucc.muts.uniq.txt", row.names = F, col.names = T, sep = "\t", quote = F)


locdiff <- fread("../10.NucleRes/sysucc.crc.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
names(locdiff) <- c("ID", "Freq", "Ratio", "Gene", "Mut", "WTscore", "Mutscore")
locdiff$DeltaScore <- locdiff$Mutscore - locdiff$WTscore


enzyme <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F, fill=TRUE)
oncotsg <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
enzyme <- enzyme[, c(5,8,1,1)]
oncotsg <- oncotsg[, c(5,8,1,1)]
tfs <- tfs[, c(5,8,1,1)]
enzyme <- rbind(enzyme, oncotsg,tfs)

locdiff <- merge(locdiff, enzyme[, c(3, 2)], by.x = "Gene", by.y = "Gene", all.x = T)


selmuts <- c("KAT8", "DEAF1", "NFE2L1", "CMAS", "NHLRC1", "PTEN", "CHFR")
sellocdiff <- locdiff[locdiff$Gene %in% selmuts, ]



###
locdiff <- fread("../10.NucleRes/sysucc.crc.selrescount.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
names(locdiff) <- c("ID", "Freq", "Ratio", "Gene", "Mut", "WTscore", "Mutscore")
locdiff$DeltaScore <- locdiff$Mutscore - locdiff$WTscore

enzyme <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Enzymes.tsv", header = T, stringsAsFactors = F, data.table = F, fill=TRUE)
oncotsg <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Cancer-related.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
enzyme <- enzyme[, c(5,8,1,1)]
oncotsg <- oncotsg[, c(5,8,1,1)]
tfs <- tfs[, c(5,8,1,1)]
enzyme <- rbind(enzyme, oncotsg,tfs)

locdiff <- merge(locdiff, enzyme[, c(3, 2)], by.x = "Gene", by.y = "Gene", all.x = T)
locdiff <- unique(locdiff)



########pathway_enrichment########
mutnuc <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
genesmerge <- fread("../4.panMut/sig.delta.merged.txt", header = T, stringsAsFactors = F, data.table = F)

#####enriched pathways
hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

ewp1 <- enricher(genesmerge[genesmerge$Class == "Nucleus", ]$gene, pvalueCutoff = 1,
                 TERM2GENE = hall, qvalueCutoff = 1, minGSSize = 1)
p1 <- dotplot(ewp1, showCategory=10) + ggtitle("")
ggsave("../10.NucleRes/mut.nucleus.pathway.pdf", p1, width = 6, height = 3)
