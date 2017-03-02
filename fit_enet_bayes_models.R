library(data.table)
library(glmnet)
library(cvTools)

options(warn=-1)
setwd('/Users/nicoleferraro/Documents/Stanford/Winter_1617/BMI217/Project/PredictionAddiction/')

cell_line_data = fread('../CCLE/CCLE_sample_info_file_2012-10-18.txt', header=T)
drug_data = fread('../CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv', header=T)
gen_data = fread('CCLE_data_format_mutation.txt', header=T)
gen_cell_data = gen_data[,3:ncol(gen_data)]
gen_ft_data = gen_data[,1:2]
drugs = unlist(unique(drug_data[,3]))
genes = unlist(unique(gen_data[,1]))
cell_lines = intersect(unlist(unique(drug_data[,1])), unlist(colnames(gen_cell_data)))
cinds = which(unlist(cell_line_data[,1]) %in% cell_lines)
pred_data_cells = cell_line_data[cinds,]
dinds = which(unlist(drug_data[,1]) %in% cell_lines)
drug_data_cells = drug_data[dinds, ]

all_r2s = c()
for (drug in drugs) {
  dinds = which(drug_data_cells[,3] == drug)
  cell_labels = unlist(drug_data_cells[dinds,1])
  act_areas = unlist(drug_data_cells[dinds,13])
  pred_data = as.data.frame(rbind(cell_labels, act_areas))
  colnames(pred_data) = cell_labels
  pred_data = pred_data[-1,]
  ginds = which(colnames(gen_cell_data) %in% cell_labels)
  gen_drug_data = as.data.frame(gen_cell_data)[,ginds]
  pred_num_data = apply(pred_data, 1, function(x) as.numeric(x))
  flds = cvFolds(ncol(gen_drug_data), K=10)
  tot_r2 = 0
  all_coefs = gen_ft_data
  cors = apply(gen_drug_data, 1, function(x) cor(as.numeric(x), pred_num_data))
  kinds = which(cors > 0.1)
  gen_drug_data_keep = gen_drug_data[kinds,]
  for (i in 1:10) {
    train_inds <- flds$subsets[flds$which != i]
    test_inds <- flds$subsets[flds$which == i]
    pred_train = pred_num_data[train_inds]
    gen_train = gen_drug_data_keep[,train_inds]
    pred_test = pred_num_data[test_inds]
    gen_test = gen_drug_data_keep[,test_inds]
    fit = cv.glmnet((t(data.matrix(gen_train))), as.vector(pred_train), alpha=0.5)
    betas = rep(0, nrow(gen_ft_data))
    coefs = coef(fit, s="lambda.min")
    betas[coefs@i] = coefs@x
    all_coefs = cbind(all_coefs, betas)
    pred = predict(fit, t(data.matrix(gen_test)), s="lambda.min")
    r = cor(pred, pred_test, method='pearson')
    tot_r2 = tot_r2 + r^2
  }
  out_file = paste(drug, '_with_mutations.txt', sep='')
  write.table(all_coefs, out_file, sep='\t', row.names=F)
  all_r2s = c(all_r2s, tot_r2/10)
}

pred_outcomes = cbind(drugs, all_r2s)
write.table(pred_outcomes, 'drug_r2_with_mutations.txt', sep='\t', row.names=F, col.names=F)

