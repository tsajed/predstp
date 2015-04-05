###########################################################################
##  Project = Knottin Protein Project                 ##
##  Contributors = Ashiqul Islam Mishu and Tanvir Sajed         ##
##  KnottinAnalysis.R -> Creates iterative sampling of 560 proteins   ##
##    of which 156 is knottin, and rest 404 are non-knottin. Using    ##
##    svm() and predict() function it tries to measure              ##
##    efficiency of prediction by training 100 proteins, and test   ##
##    the rest randomly.                        ##
##  To run this file type : source("location")                          ##                    
##########################################################################

library("class")
library("e1071")
library("ROCR")

dataset <- read.csv("BeanWrite.txt")

push <- function(vec, item) {
vec=substitute(vec)
eval.parent(parse(text = paste(vec, ' <- c(', vec, ', ', item, ')', sep = '')), n = 1)
}
count <- 0
#To calculate percentages of knottin sequences being correctly predicted
seq1 <- seq(1:1)
seq2 <- seq(1:2)
predictions <- c()
labels <- c() 
Eff <- matrix(seq2, nrow=2, ncol=2)
numberOfSVM <- matrix(seq1, nrow=1, ncol=1)
probSVM <- matrix(seq1, nrow=1, ncol=1)

for(i in 1:1) 
{
  numberOfSVM[i] <- 0
  probSVM[i] <- 0
}

for(i in 1:1) 
{  
	test.rows <- sample(nrow(dataset),100)
	train.set <- dataset[1:537,]
	test.set <- dataset[538:nrow(dataset),]
	#classifier <- naiveBayes(train.set[,1:(ncol(dataset)-1)], train.set[,ncol(dataset)], na.action=na.omit)
        
        #tuned <- tune.svm(X28~., data = train.set, gamma = 10^(-6:-1), cost = 10^(-1:1), type="C-classification")
  #Naiveclassifier <- naiveBayes(train.set[,1:ncol(dataset)-1], train.set[,ncol(dataset)], na.action=na.omit)
        #tuned <- tune.svm(X11~., data = train.set, gamma = 10^(-6:-1), cost = 10^(-1:1), type="C-classification")

        SVMclassifier <- svm(X11 ~ ., data=as.matrix(train.set), kernel="radial", gamma=0.1, cost=0.1, type="C-classification",
                             probability=TRUE)
  #df <- data.frame(predict(Naiveclassifier, test.set[,1:ncol(dataset)-1], type='raw'), test.set[,ncol(dataset)])
        
        #SVMprediction <- predict(SVMclassifier, as.matrix(test.set[,-ncol(dataset)]), probability=TRUE)
        #tab <- table(pred = SVMprediction, true = test.set[,ncol(dataset)])

        SVMPredict <- predict(SVMclassifier, as.matrix(test.set[,-ncol(dataset)]), probability=TRUE)
        prob <- attr(SVMPredict, 'probabilities')
        prediction <- as.character(SVMPredict)
        rawresult <- cbind(prediction, prob)
        

        for(i in 1:nrow(rawresult))
        {
          if(rawresult[i,1] == 1) {
            push(numberOfSVM, i)
            push(probSVM, rawresult[i,2])
            count <- count + 1
          }
        }
        #to accomodate for empty predictions , no positives, makes
        #an matrix to an array for Perl proper conversion. need to
        #check for this when dealing with empty predictions
        if(nrow(rawresult) == 0) {
          push(numberOfSVM,0)
        }
}
print(count)
