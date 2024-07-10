#setwd("/home/yukai/projects/CBAmodel/subCellLoc/bin")
args = commandArgs(trailingOnly=TRUE)
locs = args[1]
cvnum = as.numeric(args[2])
library(data.table)

filepos <- intersect(list.files("../0.data/", pattern = locs, full.names = F),
                     intersect(list.files("../0.data/", pattern = "0.2", full.names = F),
                               list.files("../0.data/", pattern = "positive_data", full.names = F)))
fileneg <- intersect(list.files("../0.data/", pattern = locs, full.names = F),
                     intersect(list.files("../0.data/", pattern = "0.2", full.names = F),
                               list.files("../0.data/", pattern = "negative_data", full.names = F)))

posdata <- fread(file.path("../0.data/", filepos), header = F, stringsAsFactors = F, data.table = F)
names(posdata) <- c("ID", "Orga", "Locs", "Seq")
posdata$Type <- "Pos"
negdata <- fread(file.path("../0.data/", fileneg), header = T, stringsAsFactors = F, data.table = F)
names(negdata) <- c("ID", "Orga", "Locs", "Seq")
negdata$Type <- "Neg"

newdata <- rbind(negdata, posdata)
newdata$Class = ifelse(newdata$Type == "Pos", 1, 0)

newdata_nuc <- newdata[newdata$Class == 1, ]
newdata_no <- newdata[newdata$Class == 0, ]
nn_nuc1 = nrow(newdata_nuc)*0.2
nn_no1 = nrow(newdata_no)*0.2

set.seed(cvnum)
sub_nuc1 = sample(1:nrow(newdata_nuc), nn_nuc1, replace=F)
valid_nuc = newdata_nuc[sub_nuc1, ]
train_nuc = newdata_nuc[-sub_nuc1, ]

sub_no1 = sample(1:nrow(newdata_no), nn_no1, replace=F)
valid_no = newdata_no[sub_no1, ]
train_no = newdata_no[-sub_no1, ]

write.table(rbind(train_nuc, train_no), paste("../0.data/1.dataset.subloc.", locs, ".txt.train.cv", cvnum, sep = ""), row.names = F, 
            col.names = T, sep = "\t", quote = F)
write.table(rbind(valid_nuc, valid_no), paste("../0.data/1.dataset.subloc.", locs, ".txt.valid.cv", cvnum, sep = ""), row.names = F, 
            col.names = T, sep = "\t", quote = F)






