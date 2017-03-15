# Assuming you have a folder with all 24 drug models
# Inputs: joint, copy_xlfile2

#extract features
extract_features <- function(drugname,jointcop){
  cosmicids <- copy_xlfile2[copy_xlfile2$Drug_Name == drugname,]$COSMIC_ID
  featureset <- jointcop[,colnames(jointcop) %in% cosmicids]
  cell_lines <- cosmicids[cosmicids %in% colnames(jointcop)]
  return(list(featureset, cell_lines))
}

#extract labels
extract_labels <- function(drug_name,cosmic_id){
  auc <- filter(copy_xlfile2, COSMIC_ID %in% cosmic_id & Drug_Name == drug_name)
  aoc <- sapply(auc$AUC, function(x) return(1-x))
  return(aoc)
}
#run
drug_list <- c("Nilotinib", "Erlotinib", "PHA-665752", "Paclitaxel", "Sorafenib", "Lapatinib")
joint_list <- list(joint_copy1, joint_copy2, joint_copy3, joint_copy4, joint_copy5, joint_copy6)
for(i in 1:length(drug_list)){
  drug <- drug_list[i]
  featAndCL <- extract_features(drug, joint_list[[i]])
  labels <- as.vector(extract_labels(drug,featAndCL[[2]]))
  pred <- predict(model_list[[i]], features, s = "lambda.min")
  labels <- labels[!is.na(pred)]
  pred <- na.omit(pred)
  rmse <- sqrt(mean((pred-labels)^2))
  r = cor(pred, labels, method ='pearson')
  print(drug)
  print(NA%in%labels)
  print(abs(r))
}
