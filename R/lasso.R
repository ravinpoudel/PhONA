#' Run Lasso Models
#'
#' This function takes in phenotype information and count data to run lasso regression. For hyper paramter optimization, 5 fold cv method is applied.
#' @param Phenotype A phenotype variable to use for lasso regression
#' @return odata A count data for feature, with rows as sample, and columns as features/ASVs/OTUs
#' @export

lasso <- function(Phenotype, odata) {

  odata_with_pheno <- as.data.frame(odata)

  odata_with_pheno["Phenotype"] <- Phenotype

  myGrid <- expand.grid(
    alpha = 1, ## lasso
    lambda = seq(0.0001, 1, length = 100)
  )

  trControl = trainControl(
    method = "cv",
    number = 5,  # changed to 5, becasue of small data set - 10 would give each fold only 2 observations, not 4.
    verboseIter = FALSE) # TRUE prints



  model <- train(Phenotype ~., data = odata_with_pheno, method = "glmnet",
                 tuneGrid = myGrid, trControl = trControl)



  aa = data.frame(varImp(model,scale=F)$importance)
  aa["otu"] <- rownames(aa)
  bb = aa[aa$Overall > 0, ]  ### remove one with zero



  cdata = odata_with_pheno[, !colnames(odata_with_pheno) %in% "Phenotype"]

  sel_cdata = cdata[, colnames(cdata) %in% bb$otu]
  sel_cdata["Phenotype"] <- odata_with_pheno$Phenotype

  sub_model <- glm(Phenotype ~ ., data = sel_cdata)
  #summary.lm(sub_model)
  #confint(sub_model)
  coff_df = data.frame(summary(sub_model)$coefficient)[-c(1), ] # also remove row with intercept
  #coff_df_sig = coff_df[coff_df$Pr...t.. <=0.05, ]
  #coff_df[trt] <- trt
  coff_df["otus"] <- rownames(coff_df)
  coff_df["relation"] <- coff_df$Estimate
  coff_df["pvalue"] <- coff_df$Pr...t..
  coff_df <- na.omit(coff_df) # remove rows with NA

  if (dim(coff_df)[1] > 0){
    coff_df
  } else {
    data.frame()
  }
}
