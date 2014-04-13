##########################################################################
##	Project = Knottin Protein Project									##
##	Contributors = Ashiqul Islam Mishu and Tanvir Sajed					##
##	KnottinAnalysis.R -> Creates iterative sampling of 200 proteins		##
##    of which 100 is knottin, and rest 100 are non-knottin. Using		##
##	  naiveBayes() and predict() function it tries to measure 			##
##    efficiency of prediction by training 50 proteins, and test		##
##	  the rest 50 randomly.												##
##  To run this file type : source("location")                          ## 										
##########################################################################

library("class")
library("e1071", lib.loc="C:/Knottine/")

dataset <- read.csv("C:/Users/Mishu/Desktop/For another Laptop/BeanWrite.txt")

push <- function(vec, item) {
vec=substitute(vec)
eval.parent(parse(text = paste(vec, ' <- c(', vec, ', ', item, ')', sep = '')), n = 1)
}

#To calculate percentages of knottin sequences being correctly predicted
seq1 <- seq(1:1)
numberOfBayes <- matrix(seq1, nrow=1, ncol=1)
numberOfSVM <- matrix(seq1, nrow=1, ncol=1)

Naivecount <- 0
SVMcount <- 0
for(i in 1:1) 
{
	numberOfBayes[i] <- 0
	numberOfSVM[i] <- 0
}
 
	test.rows <- sample(nrow(dataset),100)
	train.set <- dataset[1:560,]
        classifier <- naiveBayes(train.set[,1:(ncol(dataset)-1)], train.set[,ncol(dataset)], na.action=na.omit)
        SVMclassifier <- svm(X28 ~ ., data=train.set, kernel="radial", gamma=0.1, cost=1, type="C-classification")
        
        maintest.set <- dataset[561:nrow(dataset),]
        
for(i in 1:nrow(maintest.set)) 
{ 
	test.set <- maintest.set[i:i,]
	#classifier <- naiveBayes(train.set[,1:(ncol(dataset)-1)], train.set[,ncol(dataset)], na.action=na.omit)
        
        #tuned <- tune.svm(X7~., data = train.set, gamma = 10^(-6:-1), cost = 10^(-1:1), type="C-classification")
        #SVMclassifier <- svm(X28 ~ ., data=train.set, kernel="radial", gamma=0.1, cost=1, type="C-classification")
        
	NaivePredict <- predict(classifier, test.set[,1:(ncol(dataset)-1)], type='raw')
        
        SVMPredict <- data.frame(predict(SVMclassifier, test.set[,-ncol(dataset)]),test.set[,ncol(dataset)])
        
	#For the 100 left in the test set. test.rows 100 proteins
	for(j in 1:nrow(NaivePredict))
	{
		# df[i,1] the column for FALSE probability, df[i,2] = TRUE probability
		# df[i,3] the correct boolean for a particular protein
		# test.rows[i] contains the index for a particular protein
		# randomly sequenced in test.rows array
		
		if( NaivePredict[j,1]<0.05 && NaivePredict[j,2]>=0.95 ) {
		  Result = TRUE 
		  #print ("protein index")
		  push(numberOfBayes,i)
                  Naivecount <- Naivecount + 1
                  cat(i,file="outBayes.txt",sep="\n",append=TRUE)
		}
                if(SVMPredict[j,1] == TRUE) {
                    push(numberOfSVM,i)
                    SVMcount <- SVMcount + 1
                    cat(i,file="outSVM.txt",sep="\n",append=TRUE)
                }	
	}
}
print("Naivecount")
print(Naivecount)
print("SVMCount")
print(SVMcount)
#source("G:/Knottine/GenomeScan/1st50Volvox/KnottinAnalysis1st50.r")