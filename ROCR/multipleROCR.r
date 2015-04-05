##########################################################################
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

dataset1 <- read.csv("Feature1.txt")
dataset2 <- read.csv("Feature2.txt")
dataset3 <- read.csv("Feature3.txt")
dataset4 <- read.csv("Feature4.txt")
dataset5 <- read.csv("Feature5.txt")
dataset6 <- read.csv("Feature6.txt")


push <- function(vec, item) {
vec=substitute(vec)
eval.parent(parse(text = paste(vec, ' <- c(', vec, ', ', item, ')', sep = '')), n = 1)
}


#To calculate percentages of knottin sequences being correctly predicted
seq1 <- seq(1:1)
seq2 <- seq(1:2)
Eff <- matrix(seq2, nrow=2, ncol=2)
numberOfSVM <- matrix(seq1, nrow=1, ncol=1)
probSVM <- matrix(seq1, nrow=1, ncol=1)

for(i in 1:1) 
{
  numberOfSVM[i] <- 0
  probSVM[i] <- 0
}

Eff[1,1] <- 0
Eff[1,2] <- 0
Eff[2,1] <- 0
Eff[2,2] <- 0

plotROCR <- function(dataset, featureSet, firstplot, color) {

  predictions <- c()
  labels <- c() 

  for(i in 1:10) 
  {  
    knottin.test <- sample(1:144,50)
          nonknottin.test <- sample(145:nrow(dataset),150)
          test.rows <- c(knottin.test, nonknottin.test)
    train.set <- dataset[-test.rows,]
    test.set <- dataset[test.rows,]
          
    #Naiveclassifier <- naiveBayes(train.set[,1:ncol(dataset)-1], train.set[,ncol(dataset)], na.action=na.omit)
          #tuned <- tune.svm(X11~., data = train.set, gamma = 10^(-6:-1), cost = 10^(-1:1), type="C-classification")

          SVMclassifier <- svm(featureSet, data=train.set, kernel="radial", gamma=0.1, cost=10, type="C-classification",
                               probability=TRUE)
    #df <- data.frame(predict(Naiveclassifier, test.set[,1:ncol(dataset)-1], type='raw'), test.set[,ncol(dataset)])
          
          #SVMprediction <- predict(SVMclassifier, as.matrix(test.set[,-ncol(dataset)]), probability=TRUE)
          #tab <- table(pred = SVMprediction, true = test.set[,ncol(dataset)])

          SVMPredict <- predict(SVMclassifier, as.matrix(test.set[,-ncol(dataset)]), probability=TRUE)
          prob <- attr(SVMPredict, 'probabilities')
          prediction <- as.character(SVMPredict)
          rawresult <- cbind(prediction, prob)
          tab <- table(pred = rawresult[,1], true = as.character(test.set[,ncol(dataset)]))
          
          Eff[1,1] <- tab[1,1] + Eff[1,1]
          Eff[2,1] <- tab[2,1] + Eff[2,1]
          Eff[1,2] <- tab[1,2] + Eff[1,2]
          Eff[2,2] <- tab[2,2] + Eff[2,2]

          predictions <- c(predictions, as.numeric(rawresult[,2]))
          labels <- c(labels, test.set[,ncol(dataset)])

          for(i in 1:nrow(rawresult))
          {
            if(rawresult[i,1] == 'TRUE') {
              push(numberOfSVM, i)
              push(probSVM, rawresult[i,2])
            }
          }
          #to accomodate for empty predictions , no positives, makes
          #an matrix to an array for Perl proper conversion. need to
          #check for this when dealing with empty predictions
          if(nrow(rawresult) == 0) {
            push(numberOfSVM,0)
          }
  }

  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  if (firstplot) {
    plot(perf, col = color)
  }
  else {
    plot(perf, add = TRUE, col = color)
  }
}

#postscript("CummulativeCaptures.eps", width = 8, height = 8, 
#           horizontal = FALSE, onefile = FALSE, paper = "special", 
#          colormodel = "cmyk", family = "ComputerModern", encoding = "TeXtext.enc")

plotROCR(dataset1, X23 ~ ., TRUE, "black")
plotROCR(dataset2, X23 ~ ., FALSE, "red")
plotROCR(dataset3, X3 ~ ., FALSE, "green")
plotROCR(dataset4, X7 ~ ., FALSE, "blue")
plotROCR(dataset5, X11 ~ ., FALSE, "violet")
plotROCR(dataset6, X11 ~ ., FALSE, "cyan")

legend(0.7,0.3, 
c("FeatureSet1","FeatureSet2","FeatureSet3","FeatureSet4","FeatureSet5","FeatureSet6"),
lty=c(1,1,1,1,1,1),
lwd=c(2.5,2.5,2.5,2.5,2.5,2.5),
col=c("black","red", "green", "blue", "violet", "cyan"))



