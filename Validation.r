#setwd("/Users/DNOL4/Desktop/pset3/project/")
#import clinical data and setup to use
proj_data <- data.frame(read.csv("Cell_line_RMA_proc_basalExp.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
colnames(proj_data) = proj_data[1, ]
proj_data = proj_data[-1, ]
library(dplyr)
library(pbapply)
library(readxl)
xlfile <- read_excel("Gene_Level_CN.xlsx",sheet = 1)
xlfile2 <- read_excel("v17_fitted_dose_response.xlsx", sheet = 1)

colnames(xlfile) <- c("gene","chr","start","stop",xlfile[1,5:ncol(xlfile)])
xlfile <- xlfile[-1,]

funk <- function(x){
  # input '5,5,H,H'
  # if n1 or n2 is -1 set cnv to 0
  split <- unlist(strsplit(x,','))
  n1 <- as.numeric(split[1])
  n2 <- as.numeric(split[2])
  if(n1 == -1 & n2 == -1){
    return (0)
  }else{
    return ((n1+n2)/2)
  }

}
xlfile_copy <- xlfile
test <- pbapply(xlfile[,-c(1,2,3,4)], 2, function(x) sapply(x, funk))
xlfile_copy[,-c(1,2,3,4)] <- test
xlfile_copy[,-c(1,2,3,4)] <- lapply(xlfile_copy[,-c(1,2,3,4)],asinh)

funkstring <- function(x){
  #returns everything after "."
  split <- unlist(strsplit(x,'DATA.'))
  return (split[2])
}
colnames(proj_data)[-c(1,2)] <- lapply(colnames(proj_data)[-c(1,2)],funkstring)

#find overlapping genes
overlappingGenes <- intersect(xlfile_copy$gene, proj_data$GENE_SYMBOLS)
xlfile_copy <- filter(xlfile_copy, gene %in% overlappingGenes)
proj_data <- filter(proj_data, GENE_SYMBOLS %in% overlappingGenes)

overlappingColumns <- c("gene", "GENE_SYMBOLS", intersect(colnames(xlfile_copy), colnames(proj_data)))

xlfile_copy <- xlfile_copy[, which(names(xlfile_copy) %in% overlappingColumns)]
proj_data <- proj_data[, which(names(proj_data) %in% overlappingColumns)]
colnames(proj_data)[1] <- "gene"
joint <- rbind(proj_data, xlfile_copy)
joint <- joint[order(joint$gene),]


#read in and clean cosmic data to find cosmic ids in xlfile2

cosmic <- data.frame(read.csv("Cell_listTue Mar  7 18_12_53 2017.csv", stringsAsFactors = FALSE, sep = ",", header = FALSE))
colnames(cosmic) <- cosmic[1,]
cosmic <- cosmic[-1,]

cosmic_ids <- intersect(xlfile2$COSMIC_ID,cosmic$Name)

xlfile2$gene <- sapply(xlfile2$COSMIC_ID, function(x) cosmic[cosmic$Name == x, 2])

#find drug names and match them to drug id/cosmic id names

drugid <- read_excel("screened_compounds.xlsx", sheet = 1)

merge_drugs <- function(id){
  drugRow <- drugid[drugid$`Drug ID` == id,]
  return (drugRow$`Drug Name`)
}
xlfile2$Drug_Name <- lapply(xlfile2$DRUG_ID, merge_drugs)

copy_xlfile2 <- xlfile2

copy_xlfile2 <- filter(copy_xlfile2, Drug_Name %in% drugs)

cclexl <- read_excel("GDSC-CCLE-CTRP_conversion.xlsx")

merge_cell_lines <- function(cosmicid){
  cellline_row <- as.data.frame(cclexl)[as.numeric(cosmicid) == as.data.frame(cclexl[,1])[[1]],]
  cellline_name <- cellline_row[,3]
  return(cellline_name)
}
#"clean up" input file
joint_copy <- joint
joint_copy <- joint_copy[,c(TRUE, which(colnames(joint_copy) %in% as.data.frame(cclexl[,1])[[1]]))]
newcolnames <- c("gene", lapply(colnames(joint_copy)[-1], merge_cell_lines))
colnames(joint_copy) <- newcolnames
joint_copy <- joint_copy[,colnames(joint_copy) != "N/A"]
joint_copy <- joint_copy[,colnames(joint_copy) %in% colnames(gen_data) | colnames(joint_copy) == "gene"]
joint_copy <- joint_copy[joint_copy$gene %in% unlist(genes[kinds]),]

#read in drugs and organize
dat1 <-data.frame(read.csv("Nilotinib _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
dat2 <-data.frame(read.csv("Erlotinib _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
dat3 <-data.frame(read.csv("PHA-665752 _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
dat4 <-data.frame(read.csv("Paclitaxel _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
dat5 <-data.frame(read.csv("Sorafenib _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))
dat6 <-data.frame(read.csv("Lapatinib _betas-2.txt", stringsAsFactors = FALSE, sep = "\t", header = FALSE))

model_genes1 <- cbind(dat1$V1, dat1$V2)
model_genes2 <- cbind(dat2$V1, dat2$V2)
model_genes3 <- cbind(dat3$V1, dat3$V2)
model_genes4 <- cbind(dat4$V1, dat4$V2)
model_genes5 <- cbind(dat5$V1, dat5$V2)
model_genes6 <- cbind(dat6$V1, dat6$V2)

#clean and order to match features in training set
joint <- readRDS("joint.rds")
joint$type <- rep(c("Expression", "CNV"), nrow(joint)/2)
joint[joint$type == "Expression", -c(1, ncol(joint))] <- pbapply(joint[joint$type == "Expression", -c(1, ncol(joint))], 1, function(x) sinh(as.numeric(as.character(x))))
joint[, -c(1, ncol(joint))] <- apply(joint[, -c(1, ncol(joint))], 2, function(x) return(as.numeric(as.character(x))))
joint[, -c(1, ncol(joint))] <- pbapply(joint[, -c(1, ncol(joint))], 1, function(x) (x-mean(x))/sd(x))
subset_data <- function(data, model_genes){
  indexList <- c()
  for (row in 1:nrow(data)){
    gene <- data$gene[row]
    type <- data$type[row]
    model_row <- filter(model_genes, V1 == gene)

    if (as.character(type) %in% model_row[,2]) indexList <- c(indexList, row)
  }
  return(indexList)
}
joint_copy1 <- joint[which(joint$gene %in% model_genes1),]
joint_copy1 <- joint_copy1[subset_data(joint_copy1, dat1),]
joint_copy2 <- joint[which(joint$gene %in% model_genes2),]
joint_copy2 <- joint_copy2[subset_data(joint_copy2, dat2),]
joint_copy3 <- joint[which(joint$gene %in% model_genes3),]
joint_copy3 <- joint_copy3[subset_data(joint_copy3, dat3),]
joint_copy4 <- joint[which(joint$gene %in% model_genes4),]
joint_copy4 <- joint_copy4[subset_data(joint_copy4, dat4),]
joint_copy5 <- joint[which(joint$gene %in% model_genes5),]
joint_copy5 <- joint_copy5[subset_data(joint_copy5, dat5),]
joint_copy6 <- joint[which(joint$gene %in% model_genes6),]
joint_copy6 <- joint_copy6[subset_data(joint_copy6, dat6),]


x <- load("/Users/DNOL4/Desktop/pset3/project/validation/models.RData")
