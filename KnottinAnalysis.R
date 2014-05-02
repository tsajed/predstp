##########################################################################
##	Project = Knottin Protein Project									##
##	Contributors = Ashiqul Islam Mishu and Tanvir Sajed					##
##	KnottinAnalysis.R -> Creates iterative sampling of ## proteins		##
##    of which 156 is knottin, and rest 567 are non-knottin. Using		##
##	  svm() and predict() function it tries to measure 			        ##
##    efficiency of prediction by training 100 proteins, and test		##
##	  the rest randomly.												##
##  To run this file type : source("location")                          ## 										
##########################################################################

library("class")
library("e1071", lib.loc="C:/Knottine/")

dataset <- read.csv("C:/Users/Ashiqul/Desktop/Data Processing File/BeanWrite.txt")

#To calculate percentages of knottin sequences being correctly predicted
seq1 <- seq(1:nrow(dataset))
seq2 <- seq(1:2)
numberOfRights <- matrix(seq1, nrow=nrow(dataset), ncol=1)
numberOfWrongs <- matrix(seq1, nrow=nrow(dataset), ncol=1) 
Eff <- matrix(seq2, nrow=2, ncol=2)

Eff[1,1] <- 0
Eff[1,2] <- 0
Eff[2,1] <- 0
Eff[2,2] <- 0

for(i in 1:nrow(dataset)) 
{
	numberOfRights[i] <- 0
	numberOfWrongs[i] <- 0
}

for(i in 1:10) 
{  
	knottin.test <- sample(1:156,50)
        nonknottin.test <- sample(157:nrow(dataset),50)
        test.rows <- c(knottin.test, nonknottin.test)
	train.set <- dataset[-test.rows,]
	test.set <- dataset[test.rows,]
        
	Naiveclassifier <- naiveBayes(train.set[,1:ncol(dataset)-1], train.set[,ncol(dataset)], na.action=na.omit)
        tuned <- tune.svm(X318~., data = train.set, gamma = 10^(-6:-1), cost = 10^(-1:1), type="C-classification")

        SVMclassifier <- svm(X318 ~ ., data=train.set, kernel="radial", gamma=0.001, cost=10, type="C-classification")
	df <- data.frame(predict(Naiveclassifier, test.set[,1:ncol(dataset)-1], type='raw'), test.set[,ncol(dataset)])
        
        prediction <- predict(SVMclassifier, test.set[,-ncol(dataset)])
        tab <- table(pred = prediction, true = test.set[,ncol(dataset)])
        
        Eff[1,1] <- tab[1,1] + Eff[1,1]
        Eff[2,1] <- tab[2,1] + Eff[2,1]
        Eff[1,2] <- tab[1,2] + Eff[1,2]
        Eff[2,2] <- tab[2,2] + Eff[2,2]

	#For the 100 left in the test set. test.rows 100 proteins
	for(i in 1:100)
	{
		# df[i,1] the column for FALSE probability, df[i,2] = TRUE probability
		# df[i,3] the correct boolean for a particular protein
		# test.rows[i] contains the index for a particular protein
		# randomly sequenced in test.rows array
		
		if(df[i,2]>=0.95 && df[i,1]<=0.05 ) {
		  Result = TRUE
                  }
		else if( df[i,1]>=0.95 && df[i,2]<=0.05) {
		  Result = FALSE }
		  
		if( Result == df[i,3] ) {
		  numberOfRights[test.rows[i]] <- numberOfRights[test.rows[i]] + 1 }
		else {
		  numberOfWrongs[test.rows[i]] <- numberOfWrongs[test.rows[i]] + 1 }		
	}
}

print(Eff[1,1] * 100/(Eff[1,1] + Eff[2,1]))
print(Eff[2,2] * 100/(Eff[1,2] + Eff[2,2]))
print((Eff[2,2] + Eff[1,1]) * 100/(Eff[1,2] + Eff[2,2] + Eff[2,1] + Eff[1,1]))

TotalEfficiency <- 0
TotalKnottins <- 0
KnottinEfficiency <- 0
NonKnottinEfficiency <- 0
TotalNonKnottins <- 0

for(i in 1:nrow(dataset)) 
{
	if( numberOfRights[i] + numberOfWrongs[i] != 0)
	{
          if(i>156)
          {
            NonKnottinEfficiency <- numberOfRights[i] + NonKnottinEfficiency
            TotalNonKnottins <- (numberOfRights[i] + numberOfWrongs[i]) + TotalNonKnottins
          }
          else
          {
            KnottinEfficiency <- numberOfRights[i] + KnottinEfficiency
            TotalKnottins <- (numberOfRights[i] + numberOfWrongs[i]) + TotalKnottins
          }
	  #print("Knottin Index : ")
          #print(i)
	  #print(ProteinEfficiency)
	  #print(numberOfWrongs[i])
	}
}
  
TotalEfficiency <- (KnottinEfficiency + NonKnottinEfficiency) * 100/(TotalKnottins + TotalNonKnottins)
KnottinEfficiency <- KnottinEfficiency * 100/TotalKnottins
NonKnottinEfficiency <- NonKnottinEfficiency * 100/TotalNonKnottins
print("TotalEfficiency : ") 
print(TotalEfficiency)
print(KnottinEfficiency)
print(NonKnottinEfficiency)