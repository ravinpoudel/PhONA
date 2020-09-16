#' Run Lasso Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export

lasso <- function(x1, odata1) {

  odata_with_pheno <- as.data.frame(odata1)

  odata_with_pheno["Phenotype"] <- x1

  myGrid <- expand.grid(
    alpha = 1, ## lasso
    lambda = seq(0.0001, 1, length = 100)
  )

  trControl = trainControl(
    method = "cv",
    number = 10,
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
