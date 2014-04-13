#!/usr/bin/perl –w
use strict;
use warnings;


#AnalyzeKnottinProtein();
#ParseThroughFiles();
WriteEverySequencesToR();
#FilterGenomeSequences();

sub ParseThroughFiles {
  opendir(DIR, "C:\\Perl64\\texts");
  my @FILES = readdir(DIR);
  
  for(my $i=2; $i<@FILES; $i++) {
    #print $FILES[$i];
    AnalyzeKnottinProtein($FILES[$i]);
    print "\n++++++++++++++++++\n";
  }
}

sub FilterGenomeSequences {
  open FILE, "<", "Vcarteri_199_peptide.fas" or die $!;
  open FILE1, ">", "VCarteri_filtered.fas" or die $!;
  
  my $length = 0;
  my $sequence = "";
  my $tag = <FILE>;
  chomp($tag);
  
  while(my $line = <FILE>) {
    
    chomp($line);
    
    if(($line =~ />/) && $length < 600) {
        
      print FILE1 $tag;
      print FILE1 "\n";
      print FILE1 $sequence;
      print FILE1 "\n";
      
      $sequence = "";
      $tag = $line;
      $length = 0;
    }  
    elsif(($line =~ />/)) {
      
      $sequence = "";
      $tag = $line;
      $length = 0;
    }
    else {
      $sequence = $line.$sequence;
      $length = length($line) + $length;
    }
  }
}

sub WriteEverySequencesToR {

  open FILE, "<", "PurifiedKandNK2MSA.fas" or die $!;
  open FILE1, "<", "PurifiedKandNK2MSA.fas" or die $!;

  open FILE2, ">", "BeanWrite.txt" or die $!;
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
        
        elsif ($lineNum <= 156 ) {
          print FILE2 ",T";
        }
        
        elsif ($lineNum > 156 ) {
          print FILE2 ",F";
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