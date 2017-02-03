#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for running SVM using different combination of features

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 1) stop("Wrong number of input parameters")


path <- inputpar[1]
comb <- inputpar[2]
seeds <- inputpar[3]


library(caret)
library(dplyr)         # Used by caret
library(kernlab)       # support vector machine 
library(pROC)	       # plot the ROC curves
library(doMC)


trainSVM <- function(data, lent, featureCols=2:5, modelsNeed=c("CSC","neutral","s=1","s=2","s=3","s=5","s=10"), classCol=1, ncores=2, trainY="", subSample=FALSE, seed=1943) {
    message(paste(modelsNeed, collapse=" "))
    registerDoMC(cores = ncores)
    
    x = apply(data[,featureCols], 2, as.numeric)
    x = apply(x, 2, function(x) {
                  mean.pool = mean(x)
                  sd.pool = sd(x)
                  (x-mean.pool)/sd.pool
              })

    y = data[,classCol]
    trainX = x[(lent+1):dim(x)[1],]
    trainX = trainX[which(data$model[(lent+1):dim(data)[1]] %in% modelsNeed),]
    if (trainY == "") {
        trainY = y[(lent+1):length(y)]
        trainY = trainY[which(data$model[(lent+1):dim(data)[1]] %in% modelsNeed)]
        if ( length(modelsNeed) > 2 ) {
            trainY = sapply(trainY, function(x){if (x == "s=5" | x == "s=10" | x == "s=2" | x == "s=3")
                                                {"selection"} else {"eneutral"}})
        } else {
            trainY = sapply(trainY, function(x){ as.character(x) })
        }
        trainY = as.factor(as.character(trainY))
    } else {
        trainY = trainY
    }
    testX = trainX
    testY = trainY
    if (subSample == TRUE) {
        set.seed(seed)
        tSize = length(trainY)
        sSize = round(tSize/5)
        if (length(modelsNeed) == 2) {
            sSize = round(tSize/2)
        }
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
    cvk = 10
    if (length(modelsNeed) == 2) {
        cvk = 5
    }
    ctrl <- trainControl(method="repeatedcv",                   # 10fold cross validation
                         number=cvk,
                         repeats=5,		                # do 5 repititions of cv
                         summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                         classProbs=TRUE)

    #Train and Tune the SVM
    message("first round training")
    
    svm.tune <- train(x=trainX,
                      y= trainY,
                      method = "svmRadial",                     # Radial kernel
                      tuneLength = 9,				# 9 values of the cost function
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
        if (length(modelsNeed) > 2) {
            pred = predict.train(svm.tune, testX)
            roc = roc(testY, as.numeric(pred))
            svm.tune = list(svm.tune=svm.tune, roc = roc)
        } else {
            pred = predict.train(svm.tune, testX, type="prob")
            pred = data.frame(pred, model=testY)
            svm.tune = list(svm.tune=svm.tune, pred = pred)
        }
    }
    return(svm.tune)
}

featureComparison <- function (data, lent, combns, features, colnames, seeds, res) {    
    fsn = features[combns]
    featureCols = match(fsn, colnames)
    classCol = match("model",colnames)
    for (s in 1:length(seeds)) {
        seedn = seeds[s]
        rn = paste(fsn, collapse="_")
        rn = paste(rn, seedn, sep="_")
        message(rn)
        message(seedn)
        res[[rn]] = trainSVM(data, lent=lent, featureCols=featureCols, subSample = TRUE, seed=seedn)
    }
    return(res)
}


setwd(path)
load("stats.merged.rda")
#load("res.ica.S.rda")

features = c("fHsub","fHss","FST","KSD","rAUC")
#features = c("X1","X2")
data = stats.merged2
#data = res.ica8.S
colnames = colnames(data)
seeds = 1943:1962         #20 times each
seeds = as.numeric(seeds)
lent = 37

combns = as.numeric(strsplit(comb,"")[[1]])

featureRes = list()
featureRes = featureComparison(data, lent, combns, features, colnames, seeds, featureRes)
outfile = paste("featureICA_", comb, ".rda", sep="")
message(outfile)
save(featureRes, file=outfile)
