#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for running SVM using different combination of features

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 2) stop("Wrong number of input parameters")


path <- inputpar[1]
lent <- inputpar[2]


library(caret)
library(dplyr)         # Used by caret
library(kernlab)       # support vector machine 
library(pROC)	       # plot the ROC curves
library(doMC)

setwd(path)

load("stats.merged.rda")

features = c("fHsub","fHss","FST","KSD","rAUC")
data = stats.merged2
colnames = colnames(data)
combns2 = combn(5,2)
combns3 = combn(5,3)
combns4 = combn(5,4)
combns5 = combn(5,5)
seeds = 1943:1962      #20 times each


featureRes = list()
#featureRes = featureComparison(data, lent, combns2, features, colnames, seeds, featureRes)
#featureRes = featureComparison(data, lent, combns3, features, colnames, seeds, featureRes)
featureRes = featureComparison(data, lent, combns4, features, colnames, seeds, featureRes)
#featureRes = featureComparison(data, lent, combns5, features, colnames, seeds, featureRes)
save(res, file="featureResCombn3.rda")


featureComparison <- function (data, lent, combns, features, colnames, seeds, res) {
    for (i in 1:dim(combns)) {
        fs = as.vector(combns[,i])
        fsn = features[fs]
        featureCols = match(fsn, colnames)
        for (s in 1:length(seeds)) {
            seedn = seeds[s]
            rn = paste(fsn, collapse="_")
            rn = paste(rn,s,sep="_")
            message(rn)
            res[[rn]] = trainSVM(data, lent=lent, featureCols=featureCols, subSample = TRUE, seed=seedn)
        }
    }
    return(res)
}


trainSVM <- function(data, lent, featureCols=2:5, classCol=1, ncores=2, trainY="", subSample=FALSE, seed=1943) {
    registerDoMC(cores = ncores)
    
    x = apply(data[,featureCols], 2, as.numeric)
    x = apply(x, 2, function(x) {
                  mean.pool = mean(x)
                  sd.pool = sd(x)
                  (x-mean.pool)/sd.pool
              })

    y = data[,classCol]
    trainX = x[(lent+1):dim(x)[1],]
    if (trainY == "") {
        trainY = y[(lent+1):length(y)]
        trainY = sapply(trainY, function(x){if (x == "s=0.05" | x == "s=0.1"){"selection"} else {"eneutral"}})
        trainY = as.factor(trainY)
    } else {
        trainY = trainY
    }
    testX = trainX
    testY = trainY
    if (subSample == TRUE) {
        set.seed(seed)
        tSize = length(trainY)
        sSize = round(tSize/5)
        message(paste("testSize and trainSize:", tSize,sSize,sep=" "))
        testI = sample(tSize, sSize)
        keepI = setdiff(1:tSize, testI)
        testX = trainX[testI,]
        testY = trainY[testI]
        trainX = trainX[keepI,]
        trainY = trainY[keepI]
        message(paste("trainSize:", length(trainY), sep=" "))
    }
    
    ## SVM start
    # First pass
    set.seed(seed)
    # Setup for cross validation
    ctrl <- trainControl(method="repeatedcv",                   # 10fold cross validation
                         repeats=5,		                # do 5 repititions of cv
                         summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                         classProbs=TRUE)

    #Train and Tune the SVM
    message("first round training")
    
    svm.tune <- train(x=trainX,
                      y= trainY,
                      method = "svmRadial",                     # Radial kernel
                      tuneLength = 15,				# 9 values of the cost function
                      #tuneGrid = grid,
                      metric="ROC",
                      trControl=ctrl,
                      scaled = FALSE)

    sigma1 = as.numeric(svm.tune$bestTune["sigma"])
    message(sigma1)
    s_incre = round(sigma1/20, 2)
    message(s_incre)
    sigmaTestRange = seq(sigma1-2*s_incre, sigma1+2*s_incre, by=s_incre)
    C1 = as.numeric(svm.tune$bestTune["C"])
    message(C1)
    c1_incre = round(C1/20, 2)
    message(c1_incre)
    CTestRange = seq(C1-2*c1_incre, C1+2*c1_incre, by=c1_incre)

    # Second pass
    set.seed(seed)
    # Use the expand.grid to specify the search space	
    grid <- expand.grid(sigma = sigmaTestRange, C = CTestRange)
    
    #Train and Tune the SVM
    message("second round training: refining sigma and C")       #tend to over estimate!!!!!!!! split data for two round
    svm.tune <- train(x=trainX,
                      y= trainY,
                      method = "svmRadial",
                      #preProc = c("center","scale"),
                      metric="ROC",
                      tuneGrid = grid,
                      trControl = ctrl,
                      scaled = FALSE)
    roc = ""
    if (subSample == TRUE) {
        pred = predict.train(svm.tune, testX)
        roc = roc(testY, as.numeric(pred))
        svm.tune = list(svm.tune=svm.tune, roc = roc)
    }
    return(svm.tune)
}
