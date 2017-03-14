#AUTHOR: Nicole Ferraro
#BMI 217 Project, Winter 2016-2017
#PURPOSE: Extract values and create plots from files created by the fit_enet_bayes_models.R script
#Assumes below libraries are pre-installed
#Directory structure is complicated, needs to be changed per individual use

library(data.table)
library(ggplot2)
library(glmnet)
library(cvTools)

#Change wd as needed
setwd('/Users/nicoleferraro/Documents/Stanford/Winter_1617/BMI217/Project/PredictionAddiction/enetDrugWithMutationScaleOutput2')

#Assumes no other files in this directory besides enet model output
files = dir()

#Loop through each file and print out best predictor
for (f in files) {
  print(f)
  drug_data = fread(f, sep='\t')
  avgs = apply(drug_data[,3:ncol(drug_data)], 1, function(x) abs(mean(x)))
  mind = which(unlist(avgs) == max(unlist(avgs)))
  best_pred = drug_data[mind,1:2]
  print(best_pred)
}

#Switch wd to access file with r2 values for each drug
#Change wd as needed
setwd('/Users/nicoleferraro/Documents/Stanford/Winter_1617/BMI217/Project/PredictionAddiction/')

r2_data = fread('drug_r2_with_mutations_scale.txt', header=F)
drugs = unlist(r2_data[,1])
cor_coef = as.matrix(r2_data[,2:ncol(r2_data)])
rcor_coef = as.matrix(apply(cor_coef, c(1,2), function(x) sqrt(as.numeric(x))))
rownames(rcor_coef) = drugs
rcor_coef = rcor_coef[order(drugs),]
ses = apply(rcor_coef, 1, function(x) sd(x)/sqrt(length(x)))
means = apply(rcor_coef, 1, function(x) mean(as.numeric(x)))
plot_cors = cbind(means, ses)
rownames(plot_cors) = rownames(rcor_coef)
#Create bar plot of correlation coefficients for each drug between predictions and known values
bplot = barplot(plot_cors[,1], names.arg=unlist(rownames(plot_cors)), main="Pearson's Correlation Coefficient by Drug", ylab='Correlation Coefficient', cex.names=0.8, ylim=c(0,0.7))
par(las=2)
#Add error bars (2*standard error)
arrows(bplot,plot_cors[,1]+2*plot_cors[,2], bplot, plot_cors[,1], angle=90, code=3, length=0.05, lwd=1.5)

#Switch wd to access AUC summary file
setwd('/Users/nicoleferraro/Documents/Stanford/Winter_1617/BMI217/Project/PredictionAddiction/enetDrugWithMutationScaleOutput/')
auc_file = fread('drug_auc.txt')
drugs = auc_file[,1]
#bar plot of AUC values for each drug of plot of TP and FP rates per probability cutoff
barplot(unlist(auc_file[,2]), names.arg=unlist(drugs), main="AUC of ROC by Drug", ylab='AUC', cex.names=0.8, ylim=c(0,1))
par(las=2)

