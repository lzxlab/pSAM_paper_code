setwd("/home/yukai/projects/CBAmodel/subCellLoc/bin")

library(data.table)
library(pROC)
library(RColorBrewer)
#######################check best mdoel#######################
rocs <- fread("./nuc.all.rocs", header = F, stringsAsFactors = F, data.table = F)
rocs$ID <- paste(rocs$V1, rocs$V2, rocs$V3, rocs$V4, rocs$V5, sep = "_")
rocsmerged <- aggregate(rocs, list(rocs$ID), median)


#######################check cv aucs for all models#######################
alllocs <- fread("./locs.txt", header = F, stringsAsFactors = F, data.table = F)

for (ids in alllocs$V1) {
  color10 <- brewer.pal(10, "Spectral")
  cv01 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv10.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv02 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv11.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv03 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv12.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv04 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv13.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv05 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv14.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv06 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv15.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv07 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv16.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv08 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv17.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv09 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv18.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cv10 <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".cv19.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cvaucs <- paste(c("cv01", "cv02", "cv03", "cv04", "cv05", "cv06", "cv07", "cv08", "cv09", "cv10"),
                  c(round(auc(cv01$label, cv01$score), 4), round(auc(cv02$label, cv02$score), 4),
                    round(auc(cv03$label, cv03$score), 4), round(auc(cv04$label, cv04$score), 4),
                    round(auc(cv05$label, cv05$score), 4), round(auc(cv06$label, cv06$score), 4),
                    round(auc(cv07$label, cv07$score), 4), round(auc(cv08$label, cv08$score), 4),
                    round(auc(cv09$label, cv09$score), 4), round(auc(cv10$label, cv10$score), 4)),
                  sep = " = ")
  pdf(paste("../2.TrainValidTest/roc_data.", ids, ".cv.pdf", sep = ""), width = 4, height = 4)
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
  
}



#######################check tr-va-te aucs for all models#######################
alllocs <- fread("./locs.txt", header = F, stringsAsFactors = F, data.table = F)

for (ids in alllocs$V1){
  color10 <- c("red", "blue", "green")
  trdata <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".train.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  vadata <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".valid.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  tedata <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".test.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cvaucs <- paste(c("Train", "Valid", "Test"),
                  c(round(auc(trdata$label, trdata$score), 4), round(auc(vadata$label, vadata$score), 4),
                    round(auc(tedata$label, tedata$score), 4)),
                  sep = " = ")
  plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
  lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
  lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
  legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
  pdf(paste("../2.TrainValidTest/roc_data.", ids, ".trvate.pdf", sep = ""), width = 4, height = 4)
  plot.roc(trdata$label, trdata$score, col = color10[1], percent=TRUE)
  lines.roc(vadata$label, vadata$score, col = color10[2], percent=TRUE)
  lines.roc(tedata$label, tedata$score, col = color10[3], percent=TRUE)
  legend("bottomright", legend=cvaucs, col=color10, lwd=2, cex=1)
  dev.off()
}


#######################pair mut sampleinfos#######################
alllocs <- fread("./locs.txt", header = F, stringsAsFactors = F, data.table = F)
allpairs <- fread("../5.panpairMut/pairmut.list", header = F, stringsAsFactors = F, data.table = F)

res <- data.frame()
for (locs in alllocs$V1) {
  for (pairs in allpairs) {
    tmpres <- data.frame(ID = locs,
                         Mut = pairs)
    res <- rbind(res, tmpres)
  }
}
write.table(res, "./locs.pairmut.txt", row.names = F, col.names = F, sep = "\t", quote = F)


#######################check cutoff for all models#######################
alllocs <- fread("./locs.txt", header = F, stringsAsFactors = F, data.table = F)

for (ids in alllocs$V1){
  trdata <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".alldata.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cvaucs <- paste(c("AUC = "),
                  c(round(auc(trdata$label, trdata$score), 4)),
                  sep = " = ")
  plot.roc(trdata$label, trdata$score, percent=TRUE, print.thres = T)
  legend("bottomright", legend=cvaucs, lwd=2, cex=1)
  pdf(paste("../2.TrainValidTest/roc_data.", ids, ".all.pdf", sep = ""), width = 4, height = 4)
  plot.roc(trdata$label, trdata$score, percent=TRUE, print.thres = T)
  legend("bottomright", legend=cvaucs, lwd=2, cex=1)
  dev.off()
}

alllocs <- fread("./locs.txt", header = F, stringsAsFactors = F, data.table = F)

for (ids in alllocs$V1){
  trdata <- fread(paste("../2.TrainValidTest/roc_data.", ids, ".alldata.txt", sep = ""), header =T, stringsAsFactors = F, data.table = F)
  cvaucs <- paste(c("AUC = "),
                  c(round(auc(trdata$label, trdata$score), 4)),
                  sep = " = ")
  print(ids)
  print(trdata$score[which.max(trdata$sn-1+trdata$sp)])
}

#######################single mutation for nucleus#######################
mutfile <- fread(paste("/home/yukai/data/TCGA/mutect_tcga/TCGA-BLCA.mutect2_snv.tsv", sep = ""),
                 header = T, stringsAsFactors = F, data.table = F)
tittle <- names(mutfile)
locs <- "Nucleus"
cutoff <- 0.2
singmut <- fread(paste("../4.panMut/merged.mut.maf.score.", locs, sep = ""), header = F, stringsAsFactors = F, data.table = F)
names(singmut) <- c(tittle[2:ncol(mutfile)], "Origin_Score", "Mut_Score")
singmut$deltaScore <- singmut$Mut_Score - singmut$Origin_Score
#singmut4sig <- singmut[(singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
annovarscore <- fread("/home/yukai/data/ANNOVAR/hg38_dbnsfp30a.txt",
                      header = T, stringsAsFactors = F, data.table = F,
                      select = c("#Chr", "Start", "Ref", "Alt", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score",
                                 "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
                                 "LRT_score", "LRT_pred", "FATHMM_score", "FATHMM_pred"))
names(annovarscore)[1] <- "Chr"
annovarscore$ID <- paste(annovarscore$Chr, annovarscore$Start, annovarscore$Ref, annovarscore$Alt, sep="_")
singmut$ID <- paste(gsub("chr", "", singmut$chrom), singmut$start, singmut$ref, singmut$alt, sep="_")
singmut <- merge(singmut, annovarscore, by.x = "ID", by.y = "ID")
write.table(singmut, paste("../4.panMut/singmut.delta.", locs, sep = ""),
            row.names = F, col.names = T, sep = "\t", quote = F)

singmut <- fread(paste("../4.panMut/singmut.delta.", locs, sep = ""), header = T, stringsAsFactors = F, data.table = F)
unidescrip <- fread("/home/yukai/data/ANNOVAR/uniprot-human-function.tab",
                    header = T, stringsAsFactors = F, data.table = F)

singmut4sig <- singmut[((singmut$Origin_Score > cutoff & singmut$Mut_Score < cutoff) |
                          (singmut$Origin_Score < cutoff & singmut$Mut_Score > cutoff)) |
                         (singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
singmut4sig <- merge(singmut4sig, unidescrip, by.x = "dna_vaf", by.y = "Entry")
names(singmut4sig)[1] = "UniID"
write.table(singmut4sig, paste("../4.panMut/singmut.delta.sig.", locs, sep = ""),
            row.names = F, col.names = T, sep = "\t", quote = F)
singmut$Class <- ifelse(singmut$deltaScore < -0.1, "Sig", "Others")
singmut$FATHMM_score <- ifelse(singmut$FATHMM_score == ".", NA, singmut$FATHMM_score)
singmut$FATHMM_score <- as.numeric(singmut$FATHMM_score)
ggplot(singmut, aes(x = Class, y = FATHMM_score, color = Class))+
  geom_boxplot()+
  stat_compare_means()+
  theme_classic2()

mutafffect <- fread(paste("../4.panMut/singmut.delta.sig.", locs, sep = ""), header = T, stringsAsFactors = F, data.table = F)

kinase <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
ace <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
tfs <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
tfs <- tfs[, c(5,8,1,1)]
names(tfs) <- names(kinase)
enzyme <- rbind(kinase, ace,tfs)

mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
mutafffect <- mutafffect[, c(1,3,9,12,13,14,20,22,24,26,28,37)]
write.table(mutafffect, paste("../4.panMut/singmut.delta.sig.annote.", locs, sep = ""),
            row.names = F, col.names = T, sep = "\t", quote = F)

######test predicted hs data correlation######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
newmut <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)

prescore <- premut[, c(1, 11)]
prescore <- prescore %>% distinct(V1, .keep_all = TRUE)
newscore <- newmut[, c(1, 11)]
newscore <- newscore %>% distinct(V1, .keep_all = TRUE)

res <- merge(prescore, newscore, by.x = "V1", by.y = "V1")
corre <- cor.test(res$V11.x,res$V11.y,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p <- ggplot(res, aes(x = V11.x, y = V11.y))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle(plottitle)+
  xlab("Pre pred score")+
  ylab("New pred score")+
  theme_classic2()
p

######test predicted mut change correlation######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
premut$delta <- premut$V12 - premut$V11
premut$ID <- paste(premut$V1, premut$V7, sep = "_")
newmut <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")


prescore <- premut[, c(14, 13)]
prescore <- prescore %>% distinct(ID, .keep_all = TRUE)
newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

res <- merge(prescore, newscore, by.x = "ID", by.y = "ID")
corre <- cor.test(res$delta.x,res$delta.y,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p <- ggplot(res, aes(x = delta.x, y = delta.y))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle(plottitle)+
  xlab("Pre delta score")+
  ylab("New delta score")+
  theme_classic2()
p

######compare top delta mutations######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
premut$delta <- premut$V12 - premut$V11
premut$ID <- paste(premut$V1, premut$V7, sep = "_")
newmut <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")

prescore <- premut[, c(14, 13)]
prescore <- prescore %>% distinct(ID, .keep_all = TRUE)
newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

prescore <- prescore[order(prescore$delta), ]
newscore <- newscore[order(newscore$delta), ]

library(VennDiagram)
venn.plot <- venn.diagram(list(PreDelta=prescore$ID[1:1000],
                               NewDelta=newscore$ID[1:1000]),
                          resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),
                          fill=c("red","yellow"),
                          main="",
                          main.cex = 2, 
                          filename = NULL)

grid.draw(venn.plot)




######select prosites######
num <- 4
prosel <- data.frame(Gene = c("DEAF1", "KAT8", "NFE2L1", "CMAS", "PTEN", "CHFR", "SYNPO2"),
                     Site = c("p.K304N", "p.R144C", "p.R647C", "p.R201W", "p.R14M", "p.P255L", "p.G4E"))

newmut <- fread("../4.panMut/merged.mut.maf.score.Nucleus", header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")

newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

newscore <- newscore[newscore$ID %in% paste(prosel$Gene, prosel$Site, sep = "_"), ]
print(newscore)


newmut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")

newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

newscore <- newscore[newscore$ID %in% paste(prosel$Gene, prosel$Site, sep = "_"), ]
print(newscore)


#######################test single mutation for nucleus#######################
allnum <- c(4:10)

for (num in allnum){
  num <- 4
  mutfile <- fread(paste("/home/yukai/data/TCGA/mutect_tcga/TCGA-BLCA.mutect2_snv.tsv", sep = ""),
                   header = T, stringsAsFactors = F, data.table = F)
  tittle <- names(mutfile)
  cutoff <- 0.4
  singmut <- fread(paste("../2.Res4Nucleus/merged.mut.maf.score.Nucleus.", num, sep = ""), header = F, stringsAsFactors = F, data.table = F)
  names(singmut) <- c(tittle[2:ncol(mutfile)], "Origin_Score", "Mut_Score")
  singmut$deltaScore <- singmut$Mut_Score - singmut$Origin_Score
  #singmut4sig <- singmut[(singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
  annovarscore <- fread("/home/yukai/data/ANNOVAR/hg38_dbnsfp30a.txt",
                        header = T, stringsAsFactors = F, data.table = F,
                        select = c("#Chr", "Start", "Ref", "Alt", "SIFT_score", "SIFT_pred", "Polyphen2_HDIV_score",
                                   "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score", "Polyphen2_HVAR_pred",
                                   "LRT_score", "LRT_pred", "FATHMM_score", "FATHMM_pred"))
  names(annovarscore)[1] <- "Chr"
  annovarscore$ID <- paste(annovarscore$Chr, annovarscore$Start, annovarscore$Ref, annovarscore$Alt, sep="_")
  singmut$ID <- paste(gsub("chr", "", singmut$chrom), singmut$start, singmut$ref, singmut$alt, sep="_")
  singmut <- merge(singmut, annovarscore, by.x = "ID", by.y = "ID")
  write.table(singmut, paste("../2.Res4Nucleus/singmut.delta.Nucleus.", num, sep = ""),
              row.names = F, col.names = T, sep = "\t", quote = F)
  
  singmut <- fread(paste("../2.Res4Nucleus/singmut.delta.Nucleus.", num, sep = ""), header = T, stringsAsFactors = F, data.table = F)
  unidescrip <- fread("/home/yukai/data/ANNOVAR/uniprot-human-function.tab",
                      header = T, stringsAsFactors = F, data.table = F)
  
  singmut4sig <- singmut[((singmut$Origin_Score > cutoff & singmut$Mut_Score < cutoff) |
                            (singmut$Origin_Score < cutoff & singmut$Mut_Score > cutoff)) |
                           (singmut$deltaScore > 0.1 | singmut$deltaScore < -0.1), ]
  singmut4sig <- merge(singmut4sig, unidescrip, by.x = "dna_vaf", by.y = "Entry")
  names(singmut4sig)[1] = "UniID"
  write.table(singmut4sig, paste("../2.Res4Nucleus/singmut.delta.sig.Nucleus.", num, sep = ""),
              row.names = F, col.names = T, sep = "\t", quote = F)
  singmut$Class <- ifelse(singmut$deltaScore < -0.1, "Sig", "Others")
  singmut$FATHMM_score <- ifelse(singmut$FATHMM_score == ".", NA, singmut$FATHMM_score)
  singmut$FATHMM_score <- as.numeric(singmut$FATHMM_score)
  ggplot(singmut, aes(x = Class, y = FATHMM_score, color = Class))+
    geom_boxplot()+
    stat_compare_means()+
    theme_classic2()
  
  mutafffect <- fread(paste("../2.Res4Nucleus/singmut.delta.sig.Nucleus.", num, sep = ""), header = T, stringsAsFactors = F, data.table = F)
  
  kinase <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/Kinase.list", header = F, stringsAsFactors = F, data.table = F, fill=TRUE)
  ace <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/Acetylation.list", header = F, stringsAsFactors = F, data.table = F)
  tfs <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/protein_class_Transcription.tsv", header = T, stringsAsFactors = F, data.table = F)
  tfs <- tfs[, c(5,8,1,1)]
  names(tfs) <- names(kinase)
  enzyme <- rbind(kinase, ace,tfs)
  
  mutafffect <- merge(mutafffect, enzyme, by.x = "UniID", by.y = "V1", all.x = T)
  mutafffect <- mutafffect[, c(1,3,9,12,13,14,20,22,24,26,28,37)]
  write.table(mutafffect, paste("../2.Res4Nucleus/singmut.delta.sig.annote.Nucleus.", num, sep = ""),
              row.names = F, col.names = T, sep = "\t", quote = F)
  
}

num <- 10
######test predicted hs data correlation######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
newmut <- fread(paste("../2.Res4Nucleus/merged.mut.maf.score.Nucleus.", num, sep = ""), header = F, stringsAsFactors = F, data.table = F)

prescore <- premut[, c(1, 11)]
prescore <- prescore %>% distinct(V1, .keep_all = TRUE)
newscore <- newmut[, c(1, 11)]
newscore <- newscore %>% distinct(V1, .keep_all = TRUE)

res <- merge(prescore, newscore, by.x = "V1", by.y = "V1")
corre <- cor.test(res$V11.x,res$V11.y,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p <- ggplot(res, aes(x = V11.x, y = V11.y))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle(plottitle)+
  xlab("Pre pred score")+
  ylab("New pred score")+
  theme_classic2()
p

######test predicted mut change correlation######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
premut$delta <- premut$V12 - premut$V11
premut$ID <- paste(premut$V1, premut$V7, sep = "_")
newmut <- fread(paste("../2.Res4Nucleus/merged.mut.maf.score.Nucleus.", num, sep = ""), header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")


prescore <- premut[, c(14, 13)]
prescore <- prescore %>% distinct(ID, .keep_all = TRUE)
newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

res <- merge(prescore, newscore, by.x = "ID", by.y = "ID")
corre <- cor.test(res$delta.x,res$delta.y,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p <- ggplot(res, aes(x = delta.x, y = delta.y))+
  geom_point()+
  geom_smooth(method="lm")+
  ggtitle(plottitle)+
  xlab("Pre delta score")+
  ylab("New delta score")+
  theme_classic2()
p

######compare top delta mutations######
premut <- fread("../../ProNuclear/6.pancancer/merged.mut.maf.score", header = F, stringsAsFactors = F, data.table = F)
premut$delta <- premut$V12 - premut$V11
premut$ID <- paste(premut$V1, premut$V7, sep = "_")
newmut <- fread(paste("../2.Res4Nucleus/merged.mut.maf.score.Nucleus.", num, sep = ""), header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")

prescore <- premut[, c(14, 13)]
prescore <- prescore %>% distinct(ID, .keep_all = TRUE)
newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

prescore <- prescore[order(prescore$delta), ]
newscore <- newscore[order(newscore$delta), ]

library(VennDiagram)
venn.plot <- venn.diagram(list(PreDelta=prescore$ID[1:1000],
                               NewDelta=newscore$ID[1:1000]),
                          resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5),
                          fill=c("red","yellow"),
                          main="",
                          main.cex = 2, 
                          filename = NULL)

grid.draw(venn.plot)




######select prosites######
num <- 4
prosel <- data.frame(Gene = c("DEAF1", "KAT8", "NFE2L1", "CMAS", "PTEN", "CHFR", "SYNPO2"),
                     Site = c("p.K304N", "p.R144C", "p.R647C", "p.R201W", "p.R14M", "p.P255L", "p.G4E"))

newmut <- fread(paste("../2.Res4Nucleus/merged.mut.maf.score.Nucleus.", num, sep = ""), header = F, stringsAsFactors = F, data.table = F)
newmut$delta <- newmut$V12 - newmut$V11
newmut$ID <- paste(newmut$V1, newmut$V7, sep = "_")

newscore <- newmut[, c(14, 13)]
newscore <- newscore %>% distinct(ID, .keep_all = TRUE)

newscore <- newscore[newscore$ID %in% paste(prosel$Gene, prosel$Site, sep = "_"), ]
print(newscore)

