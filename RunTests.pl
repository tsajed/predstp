#!/usr/bin/perl
use Statistics::R;
require ServerModule;

ServerModule::TwoGenomeSequences("emptyFile","PurifiedKandNK2.fas","testtrain.fas");
my $var = ServerModule::KnottinStructureAnalysis();
if($var == 0) {
  print "The sequences you provided may contain invalid letters that are not any valid amino acids. They may contain small letter characters. Our program only analyzes capital letter valid amino acid sequences in fasta format. Thank you.";
}
else { 
  ServerModule::DataAnalysisR();
}
