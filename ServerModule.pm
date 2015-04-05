#!/usr/bin/perl

package ServerModule;

#use strict;
#use warnings;
#use Statistics::R;

my @Sequences;
my @Names;
my @Scores1;
my @Scores2;
my @Scores3;

my @PositiveBayes;
my @PositiveSVM;

my $rVar = "library('class')
library('e1071')\n";

$rVar = $rVar."dataset <- read.csv('statsForR.txt')\n";
$rVar = $rVar."push <- function(vec, item) {
vec=substitute(vec)
eval.parent(parse(text = paste(vec, ' <- c(', vec, ', ', item, ')', sep = '')), n = 1)
}\n";
$rVar = $rVar."count <- 0
#To calculate percentages of knottin sequences being correctly predicted
seq1 <- seq(1:1)
numberOfSVM <- matrix(seq1, nrow=1, ncol=1)
probSVM <- matrix(seq1, nrow=1, ncol=1)

for(i in 1:1) 
{
	numberOfSVM[i] <- 0
	probSVM[i] <- 0
}\n";

$rVar = $rVar."for(i in 1:1) 
{  
	test.rows <- sample(nrow(dataset),100)
	train.set <- dataset[1:537,]
	test.set <- dataset[538:nrow(dataset),] 
	#classifier <- naiveBayes(train.set[,1:(ncol(dataset)-1)], train.set[,ncol(dataset)], na.action=na.omit)

	SVMclassifier <- svm(X11 ~ ., data=as.matrix(train.set), kernel='radial', gamma=0.1, cost=0.1, type='C-classification', probability=TRUE)
        
	#NaivePredict <- predict(classifier, test.set[,1:(ncol(dataset)-1)], type='raw')
        
        SVMPredict <- predict(SVMclassifier, as.matrix(test.set[,-ncol(dataset)]), probability=TRUE)
        prob <- attr(SVMPredict, 'probabilities')
        prediction <- as.character(SVMPredict)
        rawresult <- cbind(prediction, prob)

	for(i in 1:nrow(rawresult))
	{
		  		                    		                		if(rawresult[i,1] == '1') {
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
}";


sub ParseThroughFiles {
  opendir (DIR,  "C:\\Users\\Tanvir\\Documents\\Knottin Project\\pdb\\non-knottins\\PDB files less than 100aa");
  my @FILES = readdir(DIR);
  
  for(my $i=0; $i<@FILES; $i++) {
    #print $FILES[$i];
    #print"\n";
    #AnalyzeKnottinProtein($FILES[$i]);
    #print "\n++++++++++++++++++\n";
  }
  shift(@FILES);
  shift(@FILES);
  return @FILES;
}

sub SequenceAlignmentPrediction() {
  
  MuscleLinedSequences("PurifiedK.fas","PredictedSequences.fas","PreAlignedSequences.fas");
  system("muscle -in PreAlignedSequences.fas -out AlignedSequences.fas");
  
  OneGenomeSequences("AlignedSequences.fas","AlignedSequences2.fas");
  open FILE, "<", "AlignedSequences2.fas" or die $!;
  
  my %CystineCount;
  my @CystineIndices;
  my @MajorIndices;
  my $number = 0;
  
  while(my $line = <FILE>) {
    if($number < 117) {
      
      if(!($line =~ />/)) {
        for(my $i=0; $i<length($line); $i++) {
          my $ith = substr($line, $i, 1);
  
          if($ith eq "C") {
            $CystineCount{$i}++;
          }
        }
        #print($CystineIndices, ",");
        $number++;
      }
      
      if($number == 117) {
        my @Keys;
        my @Values;
      
        while(my ($key,$value) = each(%CystineCount)) {
          if($value >= 60) {
            push(@MajorIndices,$key);
            print $key, "\n";
          }
        }
      }
    }
    
    else {
      my $CountMatches = 0;
      if(!($line =~ />/)) {
        for(my $i=0; $i<length($line); $i++) {
          my $ith = substr($line, $i, 1);
          
          if($ith eq "C") {
            for(my $j=0; $j<@MajorIndices; $j++) {
              if($MajorIndices[$j] == $i) {
                $CountMatches++;
              }
            }
          }
        }
        $number++;
        $CountMatches = ($CountMatches/8)*100;
        print $number, "," , $CountMatches, "\n";
      }
    }
  }
  close FILE;
}

sub DataAnalysisKnottins() {
  open FILE, ">", "PredictedSequences.fas" or die $!;
  open FILE1, "<", "BeanWrite.txt" or die $!;
  
  my $number = 1;
  my $trainset = 0;
  my @Train;
  my @Test;
  my $nb = Algorithm::NaiveBayes->new;
  my $svm = new Algorithm::SVMLight( type => 1,              # Regression model
         biased_hyperplane => 0, # Nonbiased
         kernel_type => 1 ) ;
  
  my $line = <FILE1>;
  
  while($line = <FILE1>) {
    
    if($number <= 117) {
      my @array = split ",",$line;

      $svm->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');

      
      $nb->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');
    }
    elsif($number <= 521 && $number >117){
      my @array = split ",",$line;
      
      $svm->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');

      
      $nb->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'F');
    }
    elsif ($number > 521) {
  
      my $ar = shift(@Train);
      my @array = @$ar;
      my $result = $nb->predict
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]});
      my %resulthash = %$result;
      #print($resulthash{'T'});
      if($resulthash{'T'} >= 0.95) {
        push(@PositiveBayes,$trainset);
        print FILE $Names[$trainset],"\n";
        print FILE $Sequences[$trainset],"\n";
        print $array[0],",",$array[1],",",$array[2],",",$array[3],",",$array[4],",",$array[5],"\n";
        
      }
      $trainset++;
    }
    
    if($number == 521) {
      $nb->train;
      $svm->train(@Train);
    }
    $number++;
  }
  close FILE;
  close FILE1;
}

sub DataAnalysisPerl() {
  
  open FILE, ">", "PredictedSequences.fas" or die $!;
  open FILE1, "<", "BeanWrite.txt" or die $!;
  
  my $number = 1;
  my $trainset = 0;
  my @Train;
  my @Test;
  my $nb = Algorithm::NaiveBayes->new;
  my $svm = new Algorithm::SVM(Type   => 'C-SVC',
                            Kernel => 'radial',
                            Gamma  => 64,
                            C      => 8);
  my $line = <FILE1>;
  
  while($line = <FILE1>) {
    
    if($number <= 117) {
      my @array = split ",",$line;
      
      $svm->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');

      
      $nb->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');
    }
    elsif($number <= 521 && $number >117){
      my @array = split ",",$line;
      
      $svm->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'T');
      
      $nb->add_instance
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]},  label => 'F');
    }
    else {
      my @array = split ",",$line;
      
      my $ds = new Algorithm::SVM::DataSet(Label => 'T',
                Data  => [$array[0], $array[1], $array[2], $array[3],$array[4],$array[5], $array[6]]);
      
      my $resSVM = $svm->predict
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]});
      
      my $result = $nb->predict
      (attributes => {C1 => $array[0], C2 => $array[1], C3 => $array[2], LL => $array[3],
                      B1 =>$array[4], B2 =>$array[5], B3=>$array[6]});
      my %resulthash = %$result;
      #print($resulthash{'T'});
      if($resSVM >= 0.95) {
        push(@PositiveBayes,$trainset);
        print FILE $Names[$trainset],"\n";
        print FILE $Sequences[$trainset],"\n";
        print $array[0],",",$array[1],",",$array[2],",",$array[3],",",$array[4],",",$array[5],"\n";
        
      }
      $trainset++;
    }
    
    if($number == 521) {
      $nb->train;
      $svm->train(@Train);
    }
    $number++;
  }
  close FILE;
  close FILE1;
}

sub DataAnalysisR() {
  
  open FILE, ">", "results.fas" or die $!;

  my $R = Statistics::R->new();
  $R->run($rVar);
  my $SVM = $R->get('numberOfSVM');
  my $prob = $R->get('probSVM');
    
  @PositiveSVM = @$SVM;
  shift @PositiveSVM;
   
  @ProbabilitySVM = @$prob;
  shift @ProbabilitySVM;

  #shift removes the first element from the matrices numberOfSVM, probSVM
  #They become arrays when numbers are pushed, but if there are no true
  #positives, matrices dont become arrays, Statistics::R cannot convert
  #R matrices to arrays. I pushed another 0 for empty predictions so
  #that the matrices become arrays, therefore 1 shift is not enough
  #to make these empty, they still have a 0 in them. Check for single element
  #arrays with first element of 0. Real R indices are never 0, they start
  #from 1, so 0 is a safe check number.
  
  if((length(@PositiveSVM) == 1) and ($PositiveSVM[0] == 0)) {
    print "There are probably no positive STP toxins in your sample<br>";
  }
  else {
    print "The sequences listed below are probably STP toxins from your sample. The position of the protein in your FASTA protein list is given as Index beside the probability of prediction<br><br>";
  }  

  #The first element in both arrays is 0. That's why
  #last sequence is printed. so skip the first element.
  for(my $i=0; $i<scalar(@PositiveSVM);$i++) {
   
    if($PositiveSVM[$i] != 0 && $PositiveSVM[$i] !~ /[a-zA-Z]/) { 
      print FILE $Names[$PositiveSVM[$i]-1],"\n";
      print FILE $Sequences[$PositiveSVM[$i]-1],"\n";
      print FILE $PositiveSVM[$i], " , ", $ProbabilitySVM[$i], "\n";

      print $Names[$PositiveSVM[$i]-1],"<br>";
      print $Sequences[$PositiveSVM[$i]-1],"<br>";
      print "Index = ",$PositiveSVM[$i]," , "; 
      printf "Probability = %.2f", $ProbabilitySVM[$i]; 
      print "<br><br>";
    }
    else {
      #splice @PositiveBayes, $i, 1;
      splice @PositiveSVM, $i, 1;
    }
  }
  close FILE;
}

sub PrintPositiveFromFile() {
  
  open FILE, "<", "outSVM.txt" or die $!;
  open FILE1, ">", "PredictedSequences.fas" or die $!;
  
  while(my $line = <FILE>) {
    
    chomp $line;
    
    print FILE1 $Names[$line-1],"\t", "BondScore1 : ", $Scores1[$line-1], "\t",
                "BondScore2 : ", $Scores2[$line-1], "\t","BondScore3 : ", $Scores3[$line-1], "\n";
    print FILE1 $Sequences[$line-1],"\n";
  }
}


# Algorithm to find characteristics of knottin that are unique
# like high number of Cystine knots, and very little distance
# between cystines.

sub KnottinStructureAnalysis {
  
  open FILE, "<", "testtrain.fas" or die $!;
  open FILE1, ">", "statsForR.txt" or die $!;
  
  my $tag = <FILE>;
  my $number = 0;
  my $avOneBond = 0;
  my $avTwoBond = 0;
  my $avThreeBond = 0;
  my $avRatioCystines = 0;
  my $averageDistanceCystines = 0;

  print FILE1 "0,1,2,3,4,5,6,7,8,9,10,11\n";
  
  while(my $line = <FILE>) {
    
    if($number <= 99999999) {
      
      chomp($line);
      my $ProteinLength = 0;
      my $numberOfCystines = 0;
      my $firstIndex = -1;
      my $arrayOfCystines = 0;
      my $dashes = 0;
      my @CystineLengths ;
      my @AminoAcids;
      my @BondLengths;
      my @SortedLengths;
      my %AminoCounts;
      
      $AminoCounts{'V'} = 0;
      $AminoCounts{'I'} = 0;
      $AminoCounts{'M'} = 0;
      $AminoCounts{'G'} = 0;
      $AminoCounts{'F'} = 0;
      $AminoCounts{'A'} = 0;
      $AminoCounts{'Q'} = 0;
      $AminoCounts{'L'} = 0;
      $AminoCounts{'H'} = 0;
      $AminoCounts{'R'} = 0;
      $AminoCounts{'K'} = 0;
      $AminoCounts{'T'} = 0;
      $AminoCounts{'P'} = 0;
      $AminoCounts{'E'} = 0;
      $AminoCounts{'D'} = 0;
      $AminoCounts{'N'} = 0;
      $AminoCounts{'Y'} = 0;
      $AminoCounts{'W'} = 0;
      $AminoCounts{'S'} = 0;
      $AminoCounts{'U'} = 0;
      $AminoCounts{'O'} = 0;
      $AminoCounts{'J'} = 0;

      push(@CystineLengths,0);
      
      if(!($line =~ />/)) {
        for(my $i=0; $i<length($line); $i++) {
          my $ith = substr($line, $i, 1);
  
          if($ith eq "-") {
            $dashes++;
            $ProteinLength--;
          }
          elsif($ith eq "C") {

            if($firstIndex == -1) {
              $firstIndex = $i;
              $dashes = 0;
              #$CystineLengths[$numberOfCystines] = $i;
            }
            
            else {

              $arrayOfCystines += ($i - $firstIndex - $dashes);
              push(@CystineLengths, ($i - $firstIndex - $dashes));
              my $cys = ($i - $firstIndex - $dashes);
              
              $firstIndex = $i;
              $dashes = 0;
              #$CystineLengths[$numberOfCystines] = $i;
            }
            $numberOfCystines++;
          }
          
          elsif($ith eq "V") {
            $AminoCounts{'V'}++;
          }
          elsif($ith eq "I") {
            $AminoCounts{'I'}++;
          }
          elsif($ith eq "M") {
            $AminoCounts{'M'}++;
          }
          elsif($ith eq "G") {
            $AminoCounts{'G'}++;
          }
          elsif($ith eq "Q") {
            $AminoCounts{'Q'}++;
          }
          elsif($ith eq "F") {
            $AminoCounts{'F'}++;
          }
          elsif($ith eq "L") {
            $AminoCounts{'L'}++;
          }
          elsif($ith eq "A") {
            $AminoCounts{'A'}++;
          }
           elsif($ith eq "H") {
            $AminoCounts{'H'}++;
          }
          elsif($ith eq "R") {
            $AminoCounts{'R'}++;
          }
          elsif($ith eq "K") {
            $AminoCounts{'K'}++;
          }
          elsif($ith eq "E") {
            $AminoCounts{'E'}++;
          }
          elsif($ith eq "D") {
            $AminoCounts{'D'}++;
          }
          elsif($ith eq "T") {
            $AminoCounts{'T'}++;
          }
          elsif($ith eq "Y") {
            $AminoCounts{'Y'}++;
          }
          elsif($ith eq "S") {
            $AminoCounts{'S'}++;
          }
          elsif($ith eq "N") {
            $AminoCounts{'N'}++;
          }
          elsif($ith eq "P") {
            $AminoCounts{'P'}++;
          }
          elsif($ith eq "W") {
            $AminoCounts{'W'}++;
          }
          elsif($ith eq "X") {
	   
          }
          elsif($ith eq "O") {

          }
          elsif($ith eq "J") {

          }
          elsif($ith eq "U") {

          }

          elsif($ith eq "B") {
            $AminoCounts{'N'}++;
	        }
          elsif($ith eq "Z") {
	          $AminoCounts{'Q'}++;
	        }
          # Very important that this be reverted once paper is accepted
          elsif(($ith eq "\r") or ($ith eq '') or ($ith eq "\n")) {
	        #  $ProteinLength--;
	        }
	        else {
            return 0;
          }
          
          if(!($ith eq '-') && $firstIndex != -1) {
            push(@AminoAcids,$ith);
          }
          $ProteinLength++;
        }
        
       for(my $j=scalar(@AminoAcids)-1; $j>=0; $j--) {
          if(!($AminoAcids[$j] eq 'C')) {
            pop(@AminoAcids);
          }
          else {
            last;
          }
        }
        if(scalar(@AminoAcids) == 0) {push(@AminoAcids,'-');}
        my @Classifiers = FeatureClassification(\@AminoAcids);
        my $LeastLengthRatio = 0;
        #print($CystineIndices, ",");
        $number++;
        my $av1 = 0;
        #if($numberOfCystines > 1) {
       #   $av1 = $arrayOfCystines/($numberOfCystines - 1);
        #}
        my $ratio = 0;
        #if($CystineLengths[4] != 0) {
        #  $ratio = $CystineLengths[1]/$CystineLengths[4];
        #}

        if(scalar(@CystineLengths) == 1) {
          @BondLengths = (0,0,0);
        }

        else {
          shift(@CystineLengths);
          @SortedLengths = sort{$a <=> $b}(@CystineLengths);
          $LeastLengthRatio = $SortedLengths[0]/$ProteinLength;
          #print($CystineLengths[4]) ;
          @BondLengths = BondLengthClassifier(\@SortedLengths, \@CystineLengths);
        }
        
=cut        for(my $i=0; $i<=2; $i++){
          if(($BondLengths[$i] >= 30 || $BondLengths[$i] <= 10) && $i==0) {
            $BondLengths[1] = 0;
            $BondLengths[2] = 0;
            $BondLengths[0] = 0;
          }
          elsif(($BondLengths[$i] > 25 || $BondLengths[$i] < 4) && $i==1) {
            $BondLengths[1] = 0;
            $BondLengths[2] = 0;
            $BondLengths[0] = 0;
          }
          elsif(($BondLengths[$i] >= 30 || $BondLengths[$i] < 4 ) && $i==2) {
            $BondLengths[1] = 0;
            $BondLengths[2] = 0;
            $BondLengths[0] = 0;
          }
        }
        #if($BondLengths[0] != 0) {
          #$ratio = $ProteinLength/$BondLengths[0];
          if($ProteinLength >= 150) {
            $BondLengths[1] = 0;
            $BondLengths[2] = 0;
            $BondLengths[0] = 0;
          }
        #}
        
        if($numberOfCystines > 12) {
          $BondLengths[1] = 0;
          $BondLengths[2] = 0;
          $BondLengths[0] = 0;
=cut        }

        for(my $i=0; $i<=2; $i++){
          if($BondLengths[$i] == 0) {
            $BondLengths[$i] = -10;
          }
        }

        
        my $score1 = 18.2 - $BondLengths[0] ;
        my $score2 = 15.9 - $BondLengths[1] ;
        my $score3 = 15.7 - $BondLengths[2] ;
        if($score1 < 0) {
          $score1 = -$score1;
        }
        $score1 = 100/($score1  + 10) ;
        if($score2 < 0) {
          $score2 = -$score2;
        }
        $score2 = 100/($score2 + 10);
        if($score3 < 0) {
          $score3 = -$score3;
        }
        $score3 = 100/($score3  + 10) ;
        if($number > 521) {
          push @Scores1, $score1;
          push @Scores2, $score2;
          push @Scores3, $score3;
        }
        
        if($number <= 144 ) {
           # print($BondLengths[0],',', $BondLengths[1],',', $BondLengths[2], ',' ,$number,
            #      ',',$score1,',',$score2,',',$score3);
            #print("\n");
            $avOneBond = $avOneBond + $BondLengths[0];
            $avTwoBond = $avTwoBond + $BondLengths[1];
            $avThreeBond = $avThreeBond + $BondLengths[2];
        }
        elsif($number == 145) {
            #print $avOneBond/156, "\n";
            #print $avTwoBond/156, "\n";
            #print $avThreeBond/156, "\n";
        }
        
        my $type = "T";
        if($number > 144) {
          $type = "F";
        }
        my $x = $ProteinLength/100;
        # DoubleCC boolean identifier if there are two CC adjacents in the protein
        my $DoubleCC = 0;
        if ($SortedLengths[0] == 1 && $SortedLengths[1] == 1) {
            $DoubleCC = 1;
        }
        # Boolean if CC exists in C4 - C5 or C5 - C6
        my $CC4or5 = 0;
        if ($BondLengths[3] == 1) {
          $CC4or5 = 1;
        }
        #Those who are non-knottins the first 6 cystine indices will be 0
        #This will have a character to all the non-knottins and the non-knottins
        #can be identified easily. Think about knottin indices of pivot, or other cystines.

        #print FILE1  $ThreedistanceClassifier,",", $FourdistanceClassifier, ",", $FivedistanceClassifier, ",",$RatioClassifier,",",$numberClassifier,",", $type, "\n";
        #print FILE1  $CystineLengths[0],",",$CystineLengths[4],",", $CystineLengths[4],",",$numberOfCystines,",", $type, "\n";
        print FILE1  #$Classifiers[0],",",$Classifiers[1],",",$Classifiers[2],",",
                     $LeastLengthRatio, ",", $ProteinLength,",", $DoubleCC, ",",
                     #$BondLengths[0] ,",", $BondLengths[1],",", $BondLengths[2], ",",
                     $score1 ,",", $score2,",", $score3, ",", $CC4or5, ",",
                     $numberOfCystines,",",#$type,"\n";
                     #$AminoCounts{'V'} ,",", $AminoCounts{'A'},",", $AminoCounts{'M'}, ",",$AminoCounts{'P'}, ",",
                     #$AminoCounts{'Q'} ,",", $AminoCounts{'G'},",", $AminoCounts{'F'}, ",",$AminoCounts{'R'}, ",",
                     #$AminoCounts{'K'} ,",", $AminoCounts{'Y'},",", $AminoCounts{'S'}, ",",$AminoCounts{'N'}, ",",
                     #$AminoCounts{'H'} ,",", $AminoCounts{'D'},",", $AminoCounts{'E'}, ",",$AminoCounts{'W'}, ",",
                     $AminoCounts{'S'} ,",", $AminoCounts{'K'},",", $AminoCounts{'R'}, ",",
                     $type, "\n";
      }    
    }
  }
  #print($avThreeCystines/150, "\n3distance\n");
  #print($avFourCystines/150, "\n4distance\n");
  #print($avFiveCystines/150, "\n5distance\n");
  #print($avRatioCystines/150, "\nratio\n");
  close FILE;
  close FILE1;
  return 1;
}

sub BondLengthClassifier {
  
  my $LeastBondIndex = 0;
  my $MaxBondIndex = 0;
  my @BondLengths ;

  my ($SortedLengths, $CystineLengths) = @_;
  my @SortedLengths = @$SortedLengths;
  my @CystineLengths = @$CystineLengths;
  
  for(my $i=0; $i<scalar(@SortedLengths); $i++) {
      if($SortedLengths[0] == $CystineLengths[$i]) {
        $LeastBondIndex = $i;
        $i = scalar(@SortedLengths);
      }
      if($SortedLengths[-1] == $CystineLengths[$i]) {
        $MaxBondIndex = $i;
      }
  }
  
  # Case - Proteins starting with CC or ending with CC domains are excluded since we
  # take CC domain as the pivot for a knot to be formed.
  # Case - By observation the least loop length is either 1, 2 or 3 for every knottin
  # proteins. Thereby cystine knots cannot be formed if the least loop length is not in
  # that range.
  # BondLengths[3] boolean identifier if C4 - C5 is adjacent or C5 - C6 is adjacent
  
  if($SortedLengths[0] <= 3 && $LeastBondIndex > 1 ) {
      
      # Case - C1 - C4 , C2 - C5, C3 - C6 bonds with the pivot being first C in least loop
      # Applicable only when pivot is around the end of cystines.
      
      if(($LeastBondIndex == scalar(@SortedLengths) - 2) && scalar(@CystineLengths) >= 5) {
        my $BondLength1 = $CystineLengths[$LeastBondIndex - 3] +
                          $CystineLengths[$LeastBondIndex - 2] +
                          $CystineLengths[$LeastBondIndex - 1];
        my $BondLength2 = $CystineLengths[$LeastBondIndex - 2] +
                          $CystineLengths[$LeastBondIndex - 1] +
                          $CystineLengths[$LeastBondIndex];
        my $BondLength3 = $CystineLengths[$LeastBondIndex - 1] +
                          $CystineLengths[$LeastBondIndex] +
                          $CystineLengths[$LeastBondIndex + 1];
                          
        push(@BondLengths,$BondLength1);
        push(@BondLengths,$BondLength2);
        push(@BondLengths,$BondLength3);
        # If there are CC in C4 - C5 or C5 - C6
        if (($CystineLengths[$LeastBondIndex] == 1) || ($CystineLengths[$LeastBondIndex + 1] == 1)) {
            push(@BondLengths, 1); 
        }
        else {
            push(@BondLengths, 0);
        }
      }
      
      # Case - C1 - C4 , C2 - C5, C3 - C6 bonds with the pivot being last C in least loop
      # Applicable for pivots around middle of all Cystines
      
      elsif (($LeastBondIndex < scalar(@SortedLengths) - 2) && scalar(@SortedLengths) >= 5 &&
             ($LeastBondIndex > 1)) {
        my $BondLength1 = $CystineLengths[$LeastBondIndex - 2] +
                          $CystineLengths[$LeastBondIndex - 1] +
                          $CystineLengths[$LeastBondIndex];
        my $BondLength2 = $CystineLengths[$LeastBondIndex - 1] +
                          $CystineLengths[$LeastBondIndex] +
                          $CystineLengths[$LeastBondIndex + 1];
        my $BondLength3 = $CystineLengths[$LeastBondIndex] +
                          $CystineLengths[$LeastBondIndex + 1] +
                          $CystineLengths[$LeastBondIndex + 2];
                          
        push(@BondLengths,$BondLength1);
        push(@BondLengths,$BondLength2);
        push(@BondLengths,$BondLength3);
        # If there are CC in C4 - C5 or C5 - C6
        if (($CystineLengths[$LeastBondIndex + 1] == 1) || ($CystineLengths[$LeastBondIndex + 2] == 1)) {
            push(@BondLengths, 1); 
        }
        else {
            push(@BondLengths, 0);
        }
      }
      
      # Case - Since about some knottins had only CC pivot at the end, need to find
      # two exceptional knottins with pivots at the end using C1-C4 formula
      
      elsif(scalar(@SortedLengths) >= 5 && $LeastBondIndex == scalar(@SortedLengths) - 1) {
        
        my $BondLength1 = $CystineLengths[$LeastBondIndex - 4] +
                          $CystineLengths[$LeastBondIndex - 3] +
                          $CystineLengths[$LeastBondIndex - 2] ;
        my $BondLength2 = $CystineLengths[$LeastBondIndex - 3] +
                          $CystineLengths[$LeastBondIndex - 2] +
                          $CystineLengths[$LeastBondIndex - 1];
        my $BondLength3 = $CystineLengths[$LeastBondIndex - 2] +
                          $CystineLengths[$LeastBondIndex - 1] +
                          $CystineLengths[$LeastBondIndex ];
                          
        push(@BondLengths,$BondLength1);
        push(@BondLengths,$BondLength2);
        push(@BondLengths,$BondLength3);
        # If there are CC in C4 - C5 or C5 - C6
        if ($CystineLengths[$LeastBondIndex] == 1) {
            push(@BondLengths, 1); 
        }
        else {
            push(@BondLengths, 0);
        }
      }
    
      # Cases to be solved - CC domains on end with 6 or more cysteines. Sometimes
      # knottins exist with CC domains not being pivot
      
      else {
        push(@BondLengths,0);
        push(@BondLengths,0);
        push(@BondLengths,0);
        push(@BondLengths,0);
      }
    }
     
  else {
      push(@BondLengths,0);
      push(@BondLengths,0);
      push(@BondLengths,0);
      push(@BondLengths,0);
  }
  
  #In case of ferrodoxins and iron-sulphur electron transport proteins
  #there is huge gap between the pivot and the cystine before.
  
  return @BondLengths;
}

sub FeatureClassification {
  my ($AA) = @_;
  
  my @AminoAcids = @$AA;
  my %Classify ;
  my %Classify2 ;
  
  my $Neutral = 0;
  my $Hydrophobic = 0;
  my $Hydrophilic = 0;
  
  my $Aliphatic = 0;
  my $Aromatic = 0;
  my $Polar = 0;
  my $Acidic = 0;
  my $Basic = 0;
  my $Unique = 0;
  
  $Classify{'G'} = "neutral";
  $Classify{'H'} = "neutral";
  $Classify{'S'} = "neutral";
  $Classify{'T'} = "neutral";
  $Classify{'Q'} = "neutral";
  $Classify{'R'} = "hydrophilic";
  $Classify{'K'} = "hydrophilic";
  $Classify{'N'} = "hydrophilic";
  $Classify{'D'} = "hydrophilic";
  $Classify{'E'} = "hydrophilic";
  $Classify{'P'} = "hydrophilic";
  $Classify{'B'} = "hydrophilic";
  $Classify{'Z'} = "hydrophilic";
  $Classify{'F'} = "hydrophobic";
  $Classify{'Y'} = "hydrophobic";
  $Classify{'L'} = "hydrophobic";
  $Classify{'I'} = "hydrophobic";
  $Classify{'A'} = "hydrophobic";
  $Classify{'M'} = "hydrophobic";
  $Classify{'C'} = "hydrophobic";
  $Classify{'W'} = "hydrophobic";
  $Classify{'V'} = "hydrophobic";
  $Classify{'X'} = "unknown";
  $Classify{'*'} = "unknown";
  $Classify{'-'} = "unknown";
  
  $Classify2{'G'} = "unique";
  $Classify2{'H'} = "basic";
  $Classify2{'S'} = "polar";
  $Classify2{'T'} = "polar";
  $Classify2{'Q'} = "polar";
  $Classify2{'R'} = "basic";
  $Classify2{'K'} = "basic";
  $Classify2{'N'} = "polar";
  $Classify2{'D'} = "acidic";
  $Classify2{'E'} = "acidic";
  $Classify2{'P'} = "unique";
  $Classify2{'F'} = "aromatic";
  $Classify2{'Y'} = "aromatic";
  $Classify2{'L'} = "aliphatic";
  $Classify2{'I'} = "aliphatic";
  $Classify2{'A'} = "aliphatic";
  $Classify2{'M'} = "polar";
  $Classify2{'C'} = "polar";
  $Classify2{'B'} = "polar";
  $Classify2{'Z'} = "polar";
  $Classify2{'W'} = "aromatic";
  $Classify2{'V'} = "aliphatic";
  $Classify2{'X'} = "unknown";
  $Classify2{'*'} = "unknown";
  $Classify2{'-'} = "unknown";
  
  foreach my $key (@AminoAcids) {
    if($Classify{$key} eq "neutral") {
      $Neutral++;
    }
    elsif($Classify{$key} eq "hydrophilic") {
      $Hydrophilic++;
    }
    elsif($Classify{$key} eq "hydrophobic") {
      $Hydrophobic++;
    }
    
    if($Classify2{$key} eq "aliphatic") {
      $Aliphatic++;
    }
    elsif($Classify2{$key} eq "aromatic") {
      $Aromatic++;
    }
    elsif($Classify2{$key} eq "polar") {
      $Polar++;
    }
    elsif($Classify2{$key} eq "acidic") {
      $Acidic++;
    }
    elsif($Classify2{$key} eq "basic") {
      $Basic++;
    }
    elsif($Classify2{$key} eq "unique") {
      $Unique++;
    }
    #print($key);
  }
  $Neutral = $Neutral*100/scalar(@AminoAcids);
  $Hydrophobic = $Hydrophobic*100/scalar(@AminoAcids);
  $Hydrophilic = $Hydrophilic*100/scalar(@AminoAcids);
  
  $Aliphatic = $Aliphatic*100/scalar(@AminoAcids);
  $Aromatic = $Aromatic*100/scalar(@AminoAcids);
  $Polar = $Polar*100/scalar(@AminoAcids);
  $Acidic = $Acidic*100/scalar(@AminoAcids);
  $Basic = $Basic*100/scalar(@AminoAcids);
  $Unique = $Unique*100/scalar(@AminoAcids);
  
  return ($Neutral, $Hydrophobic, $Hydrophilic, $Aliphatic, $Aromatic,
          $Polar, $Acidic, $Basic, $Unique);
}

sub MuscleLinedSequences {
  my ($M1, $M2, $W) = @_;
  
  open FILE1, "<", $M1 or die $!;
  open FILE2, "<", $M2 or die $!;
  open FILE, ">", $W or die $!;
  
  while(my $line = <FILE1>) {
    
    chomp $line;
    
    if(!($line =~ />/)) {
      $line =~ tr/-//d ;
      print FILE $line, "\n";
    }
    else {
      print FILE $line, "\n";
    }
  }
  
  while(my $line = <FILE2>) {
    
    chomp $line;
    
    print FILE $line, "\n";
  }
  close FILE1;
  close FILE2;
  close FILE;
}

sub TwoGenomeSequences {
  
  my ($M1, $M2, $W) = @_;
  
  open FILE, "<", $M1 or die $!;
  open FILE2, "<", $M2 or die $!;
  open FILE1, ">", $W or die $!;
  
  my $length = 0;
  my $sequence = "";
  my $number = 0;
  
  while(my $line = <FILE2>) {
    
    chomp($line);
    
    if(($line =~ />/)) {
    
      if( $number != 0 ) {
         #&& $exp !~ m/hypothetical/) {
        print FILE1 $sequence;
        print FILE1 "\n";
        #push(@Sequences,$sequence);
        print FILE1 $line;
        print FILE1 "\n";
        #push(@Names, $line);
      }
      else {
        print FILE1 $line;
        print FILE1 "\n";
        #push (@Names,$line);
      }
      $sequence = "";
      $length = 0;
      $number++;
    }
    else {
      $sequence = $sequence.$line;
      $length = length($line) + $length;
    }
  }
  print FILE1 $sequence;
  print FILE1 "\n";
  #push(@Sequences, $sequence);
  $number = 0;
  $sequence = ""; 
  
  while(my $line = <FILE>) {
    
    chomp($line);
    
    if(($line =~ />/)) {
    
      if( $number != 0 ) {
         #&& $exp !~ m/hypothetical/) {
        print FILE1 $sequence;
        print FILE1 "\n";
        push(@Sequences,$sequence);
        print FILE1 $line;
        print FILE1 "\n";
        push(@Names, $line);
      }
      else {
        print FILE1 $line;
        print FILE1 "\n";
        push (@Names,$line);
      }
      $sequence = "";
      $length = 0;
      $number++;
    }
    else {
      $sequence = $sequence.$line;
      $length = length($line) + $length;
    }
  }
  print FILE1 $sequence;
  print FILE1 "\n";
  push(@Sequences, $sequence);
  
  close FILE;
  close FILE1;
  close FILE2;
}

sub OneGenomeSequences {
  
  my ($M1, $W) = @_;
  
  open FILE, "<", $M1 or die $!;
  open FILE1, ">", $W or die $!;
  
  my $length = 0;
  my $sequence = "";
  my $number = 0;
  
  while(my $line = <FILE>) {
    
    chomp($line);
    
    if(($line =~ />/)) {
    
      if( $number != 0 ) {
         #&& $exp !~ m/hypothetical/) {
        print FILE1 $sequence;
        print FILE1 "\n";

        print FILE1 $line;
        print FILE1 "\n";

      }
      else {
        print FILE1 $line;
        print FILE1 "\n";

      }
      $sequence = "";
      $length = 0;
      $number++;
    }
    else {
      $sequence = $sequence.$line;
      $length = length($line) + $length;
    }
  }
  print FILE1 $sequence;
  
  close FILE;
  close FILE1;
}

sub WriteEverySequencesToR {

  open FILE, "<", "volvox1st300Fasta.fas" or die $!;
  open FILE1, "<", "volvox1st300Fasta.fas" or die $!;

  open FILE2, ">", "Volvox1st300Write.txt" or die $!;
  my $line = <FILE1>;
  my $index = 0;
  my $length = 0;
  
  
  # needs to find the column numbers for R column extraction
  # First sequence that does not have a '>' is analyzed
  # column by column and column numbers are printed to file
  # However if a new line comes before the total sequence
  # ends, the while loop goes to that new sequence line.
  # The while loop will end when it finds that the first
  # sequence has ended and a second sequence is about to
  # start, because it only needs a sequence to calculate
  # column numbers
  # $length keeps track of total sequence length as
  # new lines are added, and $index will be used in
  # for loop to print from $index to $length.

  while($line = <FILE1>) {
    
    chomp($line);
    $length = length($line) + $length;
    
    if(!($line =~ />/)) {
      
      for(my $i = $index; $i<$length; $i++) {
        print FILE2 $i;
        if($i != length($line) - 1) {
          print FILE2 ",";
        }
        else {
          print FILE2 ",";
          print FILE2 ($i + 1);
        }
      }
      $index = $index + $length;
    }
    
    else {
      print FILE2 "\n";
      last;
    }          
  }

  # Needs to print every sequence character by character
  # seperated by comma. It will also look for sequences
  # that are composed of multiple lines. Only when
  # a new line containing '>' is reached, a new line is
  # printed signalling the end of a previous sequence.
  
  my $lineNum = 0;
  
  while($line = <FILE>) {
    if($line =~ />/) {
      print FILE2 "\n";
      $lineNum++;
      next;
    }
    
    else {
      #print FILE2 "":
      chomp($line);
      for(my $i=0; $i<length($line); $i++) {
        my $ith = substr($line, $i, 1);
        
        print FILE2 $ith;
        
        if($i != length($line) - 1) {
          print FILE2 ",";
        }
        
        # Append to the end of amino acid sequence, the
        # class for every knottin protein. It could be
        # T or F, or A or S for spider. Could also
        # have multiple classes depending on toxicity
        # relative to knottin proteins
        
        elsif ($lineNum <= 3 ) {
          print FILE2 ",Agouti";
        }
        
        elsif ($lineNum > 3 && $lineNum <=6 ) {
          print FILE2 ",Algae";
        }
        
        elsif ($lineNum > 100  ) {
          print FILE2 ",Non-Toxin";
        }
      }
    }
  }
  
  close FILE;
  close FILE2;
  close FILE1;
}

sub AnalyzeKnottinProtein {
  
  my ($filename) = @_; 
  open FILE, "<", "C:\\Perl64\\texts\\".$filename;
  open FILE2, ">", "analysis.txt" or die $!;
  my @array = ();
  
  while(my $line = <FILE>) {
    if($line =~ />/) {
      next;
    }
    
    else {
      #print FILE2 "":
      chomp($line);
      my $count = FindNumberOfCystine($line);
      push @array, $count;
      
      FindDifferenceBetweenCystines($line);
      FindProlineInCystine($line);
      
      for(my $i=0; $i < length($line); $i++) {
        my $ith = substr($line, $i, 1);
        
        print FILE2 $ith;
        
        if($i != length($line) - 1) {
          print FILE2 ",";
        }
      }
      print FILE2 "\n";
    }
  }
  
  my $avCount = Average(\@array);
  print  "average ";
  print $avCount;
    
  close FILE;
  close FILE2;
}

sub FindNumberOfCystine {
  
  my($line) = @_;
  my $count = 0;

  for(my $i=0; $i < length($line); $i++) {
    my $ith = substr($line, $i, 1);
    
    if($ith eq "C") {
      $count++;
    }
  }
  print $count;
  print "\n";
  return $count;
}

sub FindDifferenceBetweenCystines {
  
  my @array = ();
  my($line) = @_;
  my $total = 0;
  my $j = 0;
  
  for(my $i=0; $i < length($line); $i++) {
    my $ith = substr($line, $i, 1);
    
    if($ith eq "C") {     
      my $difference = $i - $j;
      
      if($difference != $i) {
        push @array, $difference;
        $total = $total + $difference;
      }
      $j = $i + 1;
    }
  }
  
  #my $av = $total/@array;
  #print $av;
  #print "\n";
}

sub FindProlineInCystine {
  
  my @array = ();
  my($line) = @_;
  my $total = 0;
  my $j = 0;
  
  for(my $i=0; $i < length($line); $i++) {
    my $ith = substr($line, $i, 1);
    
    if($ith eq "C") {     
      my $difference = $i - $j;
      
      if($difference != $i) {
        push @array, $difference;
        $total = $total + $difference;
      }
      $j = $i + 1;
    }
    
    if($ith eq "P") {
      my $difference = $i - $j;
      
      if($difference != $i) {
        push @array, $difference;
        $total = $total + $difference;
      }
    }
  }
  
  #my $av = $total/@array;
  #print $av;
  #print "\n";
}

sub Average {
  
  my ($array) = @_;
  my @array = @$array;
  
  my $total = 0;
  for(my $i=0; $i<@array; $i++) {
    $total = $array[$i] + $total;
  }
  $total = $total/@array;
  return $total;
}

1;
