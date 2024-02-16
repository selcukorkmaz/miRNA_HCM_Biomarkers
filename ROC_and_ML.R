
# Library Imports
library(GEOquery)
library(limma)
library(umap)
library(dplyr)
library(caret)
library(e1071)
library(pROC)
library(patchwork)
library(cvms)

# Data Reading and Preparation

degs_gse36946_data <- read.table("data/DEGs_GSE36946.txt", header = T, sep = "\t")
dim(degs_gse36946_data)
head(degs_gse36946_data)

exprs = read.table("data/GSE36946_expression.txt", header = T, sep = "\t")

exprs2 <- cbind.data.frame(exprs[(colnames(exprs) %in% degs_gse36946_data$ID)], exprs["group"])
head(exprs2)

# ROC Curve Analysis
layout(matrix(c(1,2,0,3), 2, 2, byrow = TRUE))
roc_ILMN_3167552 <- plot(roc(exprs2$group, exprs2$ILMN_3167552), print.auc = FALSE, col = "blue")
roc_ILMN_3167552 <- plot(roc(exprs2$group, exprs2$ILMN_3168643), print.auc = FALSE, col = "darkgreen", print.auc.y = .45, add = TRUE)
roc_ILMN_3168756 <- plot(roc(exprs2$group, exprs2$ILMN_3168756), print.auc = FALSE, col = "red", print.auc.y = .4, add = TRUE)
roc_ILMN_3167463 <- plot(roc(exprs2$group, exprs2$ILMN_3167463), print.auc = FALSE, col = "black", print.auc.y = .35, add = TRUE)
roc_ILMN_3168464 <- plot(roc(exprs2$group, exprs2$ILMN_3168464), print.auc = FALSE, col = "brown", print.auc.y = .25, add = TRUE)

auc_ILMN_3167552 <- roc(exprs2$group, exprs2$ILMN_3167552)
ci_ILMN_3167552<- ci(auc_ILMN_3167552)
text(x = 0.3, y = 0.25, labels = paste0("hsa-miR-10a: ", round(auc(auc_ILMN_3167552), 3)," [", round(ci_ILMN_3167552[1], 3), "-", round(ci_ILMN_3167552[3], 3), "]"), col="blue")

auc_ILMN_3168643 <- roc(exprs2$group, exprs2$ILMN_3168643)
ci_ILMN_3168643<- ci(auc_ILMN_3168643)
text(x = 0.3, y = 0.20, labels = paste0("hsa-miR-10a*: ", round(auc(auc_ILMN_3168643), 3)," [", round(ci_ILMN_3168643[1], 3), "-", round(ci_ILMN_3168643[3], 3), "]"), col="darkgreen")

auc_ILMN_3168756 <- roc(exprs2$group, exprs2$ILMN_3168756)
ci_ILMN_3168756<- ci(auc_ILMN_3168756)
text(x = 0.3, y = 0.15, labels = paste0("hsa-miR-144*: ", round(auc(auc_ILMN_3168756), 3)," [", round(ci_ILMN_3168756[1], 3), "-", round(ci_ILMN_3168756[3], 3), "]"), col="red")

auc_ILMN_3167463 <- roc(exprs2$group, exprs2$ILMN_3167463)
ci_ILMN_3167463<- ci(auc_ILMN_3167463)
text(x = 0.3, y = 0.10, labels = paste0("hsa-miR-373: ", round(auc(auc_ILMN_3167463), 3)," [", round(ci_ILMN_3167463[1], 3), "-", round(ci_ILMN_3167463[3], 3), "]"), col="black")

auc_ILMN_3168464 <- roc(exprs2$group, exprs2$ILMN_3168464)
ci_ILMN_3168464<- ci(auc_ILMN_3168464)
text(x = 0.3, y = 0.05, labels = paste0("hsa-miR-514: ", round(auc(auc_ILMN_3168464), 3)," [", round(ci_ILMN_3168464[1], 3), "-", round(ci_ILMN_3168464[3], 3), "]"), col="brown")


roc_exprs = exprs[-ncol(exprs)]
auc_list = list()

for(i in 1:ncol(roc_exprs)){
  
  auc <- auc(roc(exprs2$group, roc_exprs[,i]))
  
  if(auc >= 0.80){
    
    auc_list[[i]] = cbind.data.frame(ID = names(roc_exprs)[i], auc = auc)
    
  }
  
  print(i)
}
  
auc_miRNA <- do.call(rbind.data.frame, auc_list)

gpl8179 = getGEO("GPL8179")
annotation = gpl8179@dataTable@table
probe_ids <- auc_miRNA$ID

miRNAs = annotation[annotation$SYMBOL %in% auc_miRNA$ID, c("SYMBOL", "miRNA_ID")]
colnames(miRNAs)[1] = "ID"
auc_miRNA = left_join(auc_miRNA, miRNAs, by="ID")

head(auc_miRNA)


# Machine Learning

topgenes = read.table("data/TopGenes_GSE36946.txt", header=T, sep="\t")
head(topgenes)

top_miRNAs_auc = topgenes[topgenes$ID %in%auc_miRNA$ID,]

ml_data = exprs2
ml_data2 = cbind.data.frame(exprs[,top_miRNAs_auc$ID], exprs["group"])
dim(ml_data2)
head(ml_data)
ml_data2[1:6,25:27]

head(ml_data)
dim(ml_data)


### Model 1 ####
# Convert group to a factor
ml_data$group <- as.factor(ml_data$group)

# Split the data into training and testing sets
set.seed(123) # For reproducibility
indexes <- createDataPartition(ml_data$group, p = 0.7, list = FALSE)
train_data <- ml_data[indexes, ]
test_data <- ml_data[-indexes, ]

# Train the SVM model with RBF kernel and 5-fold cross-validation for parameter optimization
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

###### SVM #########
svm_model <- train(group ~ ., data = train_data, method = "svmRadial",
                   trControl = train_control, preProcess = c("center", "scale"),
                   metric = "ROC")

# Predict on test data
predictions_svm <- predict(svm_model, test_data[-ncol(test_data)], type = "prob")
predictions_svm_class <- predict(svm_model, test_data[-ncol(test_data)], type = "raw")


###### RF #########
rf_model <- train(group ~ ., data = train_data, method = "rf",
                   trControl = train_control, preProcess = c("center", "scale"),
                   metric = "ROC")

# Predict on test data
predictions_rf <- predict(rf_model, test_data[-ncol(test_data)], type = "prob")
predictions_rf_class <- predict(rf_model, test_data[-ncol(test_data)], type = "raw")


# Calculate AUC
roc_result_svm <- plot(roc(response = test_data$group, predictor = predictions_svm[,2]), col = "blue")
roc_result_rf <- plot(roc(response = test_data$group, predictor = predictions_rf[,2]), col = "red", print.auc.y = .4, add = TRUE, lty=5)

auc_SVM <- roc(response = test_data$group, predictor = predictions_svm[,2])
ci_SVM <- ci(auc_SVM)
text(x = 0.2, y = 0.25, labels = paste0("SVM: ", round(auc(auc_SVM), 3)," [", round(ci_SVM[1], 3), "-", round(ci_SVM[3], 3), "]"), col="blue")

auc_RF <- roc(response = test_data$group, predictor = predictions_rf[,2])
ci_RF <- ci(auc_RF)
text(x = 0.2, y = 0.20, labels = paste0("RF: ", round(auc(auc_RF), 3)," [", round(ci_RF[1], 3), "-", round(ci_RF[3], 3), "]"), col="red")

  
svm_matrix = confusion_matrix(test_data[,ncol(test_data)], predictions_svm_class)  

confusionMatrix(table(predictions_svm_class, test_data[,ncol(test_data)]), positive = "HCM")  

rf_matrix = confusion_matrix(test_data[,ncol(test_data)], predictions_rf_class)  

confusionMatrix(table(predictions_rf_class, test_data[,ncol(test_data)]), positive = "HCM")  



p1 = plot_confusion_matrix(svm_matrix$`Confusion Matrix`[[1]],
                      rm_zero_percentages = FALSE,
                      rm_zero_text = FALSE) 

p2 = plot_confusion_matrix(rf_matrix$`Confusion Matrix`[[1]],
                           rm_zero_percentages = FALSE,
                           rm_zero_text = FALSE) 


### Model 2 ####
# Convert group to a factor
ml_data2$group <- as.factor(ml_data2$group)

# Split the data into training and testing sets
set.seed(123) # For reproducibility
indexes <- createDataPartition(ml_data2$group, p = 0.7, list = FALSE)
train_data <- ml_data2[indexes, ]
test_data <- ml_data2[-indexes, ]

# Train the SVM model with RBF kernel and 5-fold cross-validation for parameter optimization
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

###### SVM #########
svm_model <- train(group ~ ., data = train_data, method = "svmRadial",
                   trControl = train_control, preProcess = c("center", "scale"),
                   metric = "ROC")

# Predict on test data
predictions_svm <- predict(svm_model, test_data[-ncol(test_data)], type = "prob")
predictions_svm_class <- predict(svm_model, test_data[-ncol(test_data)], type = "raw")


###### RF #########
rf_model <- train(group ~ ., data = train_data, method = "rf",
                  trControl = train_control, preProcess = c("center", "scale"),
                  metric = "ROC")

# Predict on test data
predictions_rf <- predict(rf_model, test_data[-ncol(test_data)], type = "prob")
predictions_rf_class <- predict(rf_model, test_data[-ncol(test_data)], type = "raw")

# Calculate AUC


roc_result_svm <- plot(roc(response = test_data$group, predictor = predictions_svm[,2]), col = "blue")
roc_result_rf <- plot(roc(response = test_data$group, predictor = predictions_rf[,2]), col = "red", print.auc.y = .4, add = TRUE, lty=5)

auc_SVM <- roc(response = test_data$group, predictor = predictions_svm[,2])
ci_SVM <- ci(auc_SVM)
text(x = 0.2, y = 0.25, labels = paste0("SVM: ", round(auc(auc_SVM), 3)," [", round(ci_SVM[1], 3), "-", round(ci_SVM[3], 3), "]"), col="blue")

auc_RF <- roc(response = test_data$group, predictor = predictions_rf[,2])
ci_RF <- ci(auc_RF)
text(x = 0.2, y = 0.20, labels = paste0("RF: ", round(auc(auc_RF), 3)," [", round(ci_RF[1], 3), "-", round(ci_RF[3], 3), "]"), col="red")


svm_matrix = confusion_matrix(test_data[,ncol(test_data)], predictions_svm_class)  
rf_matrix = confusion_matrix(test_data[,ncol(test_data)], predictions_rf_class)  

confusionMatrix(table(predictions_svm_class, test_data[,ncol(test_data)]), positive = "HCM")  
confusionMatrix(table(predictions_rf_class, test_data[,ncol(test_data)]), positive = "HCM")  

p3 = plot_confusion_matrix(svm_matrix$`Confusion Matrix`[[1]],
                           rm_zero_percentages = FALSE,
                           rm_zero_text = FALSE) 

p4 = plot_confusion_matrix(rf_matrix$`Confusion Matrix`[[1]],
                           rm_zero_percentages = FALSE,
                           rm_zero_text = FALSE) 



p = (p1|p2) / (p3|p4)
p = p+ plot_annotation(tag_levels = 'A')

ggsave('~/Documents/Studies/LVH_DeepLearning/Figures/confusion_matrix.png', p, width = 12, height = 8)


