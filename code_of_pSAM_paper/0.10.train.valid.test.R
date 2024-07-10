#setwd("/home/yukai/projects/CBAmodel/subCellLoc/bin")
args = commandArgs(trailingOnly=TRUE)
locs = args[1]
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
nn_nuc2 = nrow(newdata_nuc)*0.1
nn_no1 = nrow(newdata_no)*0.2
nn_no2 = nrow(newdata_no)*0.1

set.seed(123456)
sub_nuc1 = sample(1:nrow(newdata_nuc), nn_nuc1, replace=F)
valid_nuc = newdata_nuc[sub_nuc1, ]
res_nuc = newdata_nuc[-sub_nuc1, ]
sub_nuc2 = sample(1:nrow(res_nuc), nn_nuc2, replace=F)
test_nuc = res_nuc[sub_nuc2, ]
train_nuc = res_nuc[-sub_nuc2, ]

sub_no1 = sample(1:nrow(newdata_no), nn_no1, replace=F)
valid_no = newdata_no[sub_no1, ]
res_no = newdata_no[-sub_no1, ]
sub_no2 = sample(1:nrow(res_no), nn_no2, replace=F)
test_no = res_no[sub_no2, ]
train_no = res_no[-sub_no2, ]

write.table(rbind(train_nuc, train_no), paste("../0.data/1.dataset.subloc.", locs, ".txt.train", sep = ""), row.names = F, 
            col.names = T, sep = "\t", quote = F)
write.table(rbind(valid_nuc, valid_no), paste("../0.data/1.dataset.subloc.", locs, ".txt.valid", sep = ""), row.names = F, 
            col.names = T, sep = "\t", quote = F)
write.table(rbind(test_nuc, test_no), paste("../0.data/1.dataset.subloc.", locs, ".txt.test", sep = ""), row.names = F, 
            col.names = T, sep = "\t", quote = F)






