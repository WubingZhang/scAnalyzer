trainAnnotator <- function(refobj, markers, method = "svm", nthread = 4){
  cl <- makePSOCKcluster(nthread)
  registerDoParallel(cl)
  trainDat <- FetchData(refobj, unlist(markers), slot = "scale.data")
  #### random forest will take more time (more core is recommended)
  message(Sys.time(), " Train the model ...")
  fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs = TRUE)
  if("svm" %in% tolower(method)){
    model <- train(x=trainDat, y=Idents(refobj), method = "svmLinear",
                 trControl = fitControl, tuneLength = 10)
  }else{
    model <- train(x=trainDat, y=Idents(refobj), method = "glmnet", family = "binomial",
                      trControl = fitControl, tuneLength = 10)
  }
  stopCluster(cl)
  return(model)
}
