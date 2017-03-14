#AUTHOR: Nicole Ferraro
#BMI 217 Project, Winter 2016-2017
#PURPOSE: Fit an elastic net and naive bayes model to data from CCLE and cross-validate
#Must have below libraries pre-installed

library(data.table)
library(e1071)
library(glmnet)
library(cvTools)
library(ROCR)

#Set personal working directory
setwd('/Users/nicoleferraro/Documents/Stanford/Winter_1617/BMI217/Project/PredictionAddiction/')
#Read in genes from validation set (extracted from TCGA data elsewhere)
gene_list = readRDS('TCGA-Gene-Set')

#Assumes data files are stored in CCLE folder, one directory up, change accordingly
cell_line_data = fread('../CCLE/CCLE_sample_info_file_2012-10-18.txt', header=T)
drug_data = fread('../CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv', header=T)
#Feature data created from the arrange_data.py script, change file if using different feature input set
gen_data = fread('CCLE_data_format_mutation.txt', header=T)
ccle_genes = unlist(unique(gen_data[,1]))
#Subset to genes existing in both training and test sets for validation
overlap_genes = intersect(gene_list, ccle_genes)
oinds = which(ccle_genes %in% overlap_genes)
gen_data = gen_data[oinds,]
#set if using validation set with only expression and CNV data
validate = T
if (validate == T) {
  kinds = which(gen_data[,'Type'] == 'Expression' | gen_data[, 'Type'] == 'CNV')
  gen_data = gen_data[kinds,]
}

#Split data and feature names
gen_cell_data = gen_data[,3:ncol(gen_data)]
gen_ft_data = gen_data[,1:2]
drugs = unlist(unique(drug_data[,3]))
genes = unlist(unique(gen_data[,1]))

#Extract cell lines with data in both drug files and feature files
cell_lines = intersect(unlist(unique(drug_data[,1])), unlist(colnames(gen_cell_data)))
cinds = which(unlist(cell_line_data[,1]) %in% cell_lines)
pred_data_cells = cell_line_data[cinds,]
dinds = which(unlist(drug_data[,1]) %in% cell_lines)
drug_data_cells = drug_data[dinds, ]

#Fit an elastic net regression model for each drug
#Use lambda with smallest MSE, determined via cross-validation
#Conduct 10-fold cross validation to obtain R^2 for the predicted drug activity area vs. known
#All files are written to current working directory
all_r2s = matrix(c(rep(0,10*24)), nrow=24, ncol=10)
rownames(all_r2s) = drugs
#Train using cross-validation, or train on entire dataset (set to F if validating on external dataset)
cross = F
#Only look at subset of drugs if using external test set
overlap = F
if (overlap == T) { drugs = c("Nilotinib", "Erlotinib", "PHA-665752", "Paclitaxel", "Sorafenib", "Lapatinib") }
model_list = list()
for (drug in drugs) {
  dinds = which(drug_data_cells[,3] == drug)
  cell_labels = unlist(drug_data_cells[dinds,1])
  act_areas = unlist(drug_data_cells[dinds,13])
  pred_data = as.data.frame(rbind(cell_labels, act_areas))
  colnames(pred_data) = cell_labels
  pred_data = pred_data[-1,]
  ginds = which(colnames(gen_cell_data) %in% cell_labels)
  gen_data = as.data.frame(gen_cell_data)[,ginds]
  #Normalize data via z-score
  gen_drug_data = scale(gen_data, center=TRUE, scale=TRUE)
  pred_num_data = apply(pred_data, 1, function(x) as.numeric(x))
  #For computational efficiency, only include features with r > 0.1 with outcome
  cors = apply(gen_drug_data, 1, function(x) cor(as.numeric(x), pred_num_data))
  kinds = which(cors > 0.1)
  gen_drug_data_keep = gen_drug_data[kinds,]
  gen_ft_data_keep = gen_ft_data[kinds,]
  if (cross == F) {
    fit = cv.glmnet((t(data.matrix(gen_drug_data_keep))), as.vector(pred_num_data), alpha=0.5)
    betas = rep(0, nrow(gen_ft_data_keep))
    coefs = coef(fit, s="lambda.min")
    betas[coefs@i] = coefs@x
    fit_out = cbind(gen_ft_data_keep, betas)
    out_file = paste(drug, '_betas.txt')
    write.table(fit_out, file=out_file, sep='\t', row.names=F, col.names=F)
    #Save models for validation step
    model_list = append(model_list, list(fit))
  }
  else {
    flds = cvFolds(ncol(gen_drug_data_keep), K=10)
    all_coefs = gen_ft_data_keep
    #Ten iterations of ten fold cross-validation
    for (j in 1:10) {
      tot_r2 = 0
      for (i in 1:10) {
        train_inds <- flds$subsets[flds$which != i]
        test_inds <- flds$subsets[flds$which == i]
        pred_train = pred_num_data[train_inds]
        gen_train = gen_drug_data_keep[,train_inds]
        pred_test = pred_num_data[test_inds]
        gen_test = gen_drug_data_keep[,test_inds]
        fit = cv.glmnet((t(data.matrix(gen_train))), as.vector(pred_train), alpha=0.5)
        betas = rep(0, nrow(gen_ft_data_keep))
        coefs = coef(fit, s="lambda.min")
        betas[coefs@i] = coefs@x
        all_coefs = cbind(all_coefs, betas)
        pred = predict(fit, t(data.matrix(gen_test)), s="lambda.min")
        r = cor(pred, pred_test, method='pearson')
        r2_val = r^2
        if (is.na(r2_val)) { r2_val = 0 }
        tot_r2 = tot_r2 + r2_val
      }
      all_r2s[drug, j] = tot_r2/10
    }
    out_file = paste(drug, '_with_mutations_scale.txt', sep='')
    write.table(all_coefs, out_file, sep='\t', row.names=F)
  }
}

pred_outcomes = cbind(drugs, all_r2s)
write.table(pred_outcomes, 'drug_r2_with_mutations_scale.txt', sep='\t', row.names=F, col.names=F)

#Fit a naive bayes classifier for each drug
#Thresholding is based on the median sensitivity value
for (drug in drugs) {
  dinds = which(drug_data_cells[,3] == drug)
  cell_labels = unlist(drug_data_cells[dinds,1])
  act_areas = unlist(drug_data_cells[dinds,13])
  pred_data = as.data.frame(rbind(cell_labels, act_areas))
  colnames(pred_data) = cell_labels
  pred_data = pred_data[-1,]
  #Threshold sensitivity based on the median value
  med = median(as.numeric((as.matrix(pred_data))))
  ginds = which(as.matrix(pred_data) >= med)
  linds = which(as.matrix(pred_data) < med)
  pred_data[ginds] = 1
  pred_data[linds] = 0
  pred_data = as.numeric(as.character(unlist(pred_data)))
  ginds = which(colnames(gen_cell_data) %in% cell_labels)
  gen_data = as.data.frame(gen_cell_data)[,ginds]
  #Normalize data via z-score
  gen_drug_data = scale(gen_data, center=TRUE, scale=TRUE)
  if (cross == F) {
    fit = naiveBayes((t(data.matrix(gen_drug_data_keep))), as.vector(pred_data))
    pred = predict(fit, t(data.matrix(gen_test)), type='raw')
    pr <- prediction(pred, pred_test)
    prf <- performance(pr, "tpr", "fpr")
    #Uncomment below line to see ROC plot for each drug
    #plot(prf, main=drug, cex.lab=1.4, cex.main=2)
    auc.perf = performance(pr, measure = 'auc')
    write.table(c(drug, auc.perf@y.values), file='drug_auc_testset.txt', sep='\t', row.names=F, col.names=F, append=T)
  }
  else {
    flds = cvFolds(ncol(gen_drug_data), K=10)
    #For computational efficiency, only include features with r > 0.1 with outcome
    cors = apply(gen_drug_data, 1, function(x) cor(as.numeric(x), pred_data))
    kinds = which(cors > 0.1)
    gen_drug_data_keep = gen_drug_data[kinds,]
    overall_pred = c()
    overall_pred_test = c()
    for (i in 1:10) {
      train_inds <- flds$subsets[flds$which != i]
      test_inds <- flds$subsets[flds$which == i]
      pred_train = pred_data[train_inds]
      gen_train = gen_drug_data_keep[,train_inds]
      pred_test = pred_data[test_inds]
      gen_test = gen_drug_data_keep[,test_inds]
      fit = naiveBayes((t(data.matrix(gen_train))), as.vector(pred_train))
      pred = predict(fit, t(data.matrix(gen_test)), type='raw')
      overall_pred = c(overall_pred, pred[,2])
      overall_pred_test = c(overall_pred_test, pred_test)
    }
    pr <- prediction(overall_pred, overall_pred_test)
    prf <- performance(pr, "tpr", "fpr")
    #Uncomment below line to see ROC plot for each drug
    #plot(prf, main=drug, cex.lab=1.4, cex.main=2)
    auc.perf = performance(pr, measure = 'auc')
    write.table(c(drug, auc.perf@y.values), file='drug_auc.txt', sep='\t', row.names=F, col.names=F, append=T)
  }
}
