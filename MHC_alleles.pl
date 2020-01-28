#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Slurp;#USE:
use Data::Dumper;
use List::Util qw(max);

use Cwd;#USE

my ($var_nam, $var_rev, $var_upd)= ("MHC_Alleles", "1.0.0", "20151022");
my ($raw_Data) = (""); 
my $workingPath= getcwd();
GetOptions("in|i=s"=> \$raw_Data);
my $inPath = <$workingPath/$raw_Data>;
open my $inFile, '<', "$inPath" or die "Can't open $inPath: $! ";
my $var_path = getcwd();
require"$var_path/_script/basic_script.pl";

mhc_A_B(Path =>$workingPath,
		InFile	=>$inFile);

sub mhc_A_B{
		my (%args)	=(@_);
		my $var_pat	= $args{Path}; 
		my $inFile	= $args{InFile};
		directoryStructure ( Path => $var_pat, 
			Folders		 => ["_reports"]);
		my @allSampleArray;
		my @allAlleleArray;
		my @fileArray;
		my %haploHash 	= filesToHash( PathIn => "$var_path/_resources/haploData.csv" );
		my %haploHashB	= filesToHash( PathIn => "$var_path/_resources/haploDataB.csv" );
		my @A_file;
		my @B_file;
		my @B2_file;
		my @data_file;

#********************************* Sample Number *******************************
		foreach my $inFileLine(<$inFile>){
				push @fileArray, $inFileLine;
				my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $inFileLine;
				push @allSampleArray, $sample;
				push @allAlleleArray, $allele;
		}
		shift @allSampleArray; #remove header
		shift @allAlleleArray; #remove header
		my @uniqueSampleArray = uniqueArray(@allSampleArray);
		my @uniqueAlleleArray = uniqueArray(@allAlleleArray);
		my $uniqueLength = @uniqueSampleArray;
		push @A_file,  "Total Number of Samples: ".$uniqueLength."\n";
		push @B_file,  "Total Number of Samples: ".$uniqueLength."\n";
		push @B2_file, "Total Number of Samples: ".$uniqueLength."\n";
#************************************* Iterate through each Sample ****************
		foreach my $currentSample (@uniqueSampleArray){
				push @A_file, "\n\nCurrent Sample: ".$currentSample."\n";
				push @B_file, "\n\nCurrent Sample: ".$currentSample."\n";
				push @B2_file, "\n\nCurrent Sample: ".$currentSample."\n";
				my @sampleArray;
				foreach my $line(@fileArray){
						my($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $line;
						if ($sample eq $currentSample){ push @sampleArray, $line; }
				}
		#********************* A1Array **********************************
				my @A1Array = grep { $_ =~ /A1\_/ } @sampleArray;#INPUT DATA
				if(@A1Array){
						my @revSortedArray = map{$_->[1]} sort{$b->[0]<=>$a->[0]} map{[substr($_,0,index($_,"\t")),$_]}@A1Array;
						push @A_file, "Total Allele Calls: ".@revSortedArray."\n";
						my @usedTerms;
						foreach my $element (@revSortedArray){
								my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $element;
								push @A_file, $allele."\t".$readCount."\n";
								# print "\n".$sample." ";
								my @haploDataArray = mhc_allelesToSearch(Term => "A1",
																			Allele => $allele);
								@haploDataArray	= uniqueArray(@haploDataArray);
								$allele =~ s/\&\_/\|/g;
								foreach my $term(@haploDataArray){
												# print "\n\tTerm: ".$term;
												my $noRepeat = 0;
												foreach my $used(@usedTerms){
														if($used eq $term){$noRepeat = 1;}#print "\tUsed";}
															else{}
														# print "\t".$used;
												}
												if($noRepeat == 0){
														push @usedTerms, $term;
														my @presentTermArray;
														foreach my $key (keys %{$haploHash{haploData}}){ #HAPLODATA
																if($haploHash{haploData}{$key}{'Diagnostic Major'} eq $term){
																		push @presentTermArray, $key;
																}#if
														}#key
														my $arrayLength = @presentTermArray;
														my @addAllSequences;
														if($arrayLength >1){
																foreach my $key(@presentTermArray){
																				my $increment = 0;
																				my @majorAlleles = mhc_get_majorAlleles(fileIn=> "haploData",
																															KeyIn => $key,
																															HashIn=> \%haploHash);
																				push @A_file, "\n";
																				my $A_length = @majorAlleles;
																				if($A_length >= 1){
																						foreach my $allele(@majorAlleles){
																								$allele =~ s/\*/\_/g;
																								my @assMajAllele = grep{$_ =~ /$allele/} @sampleArray;
																								if(@assMajAllele){$increment++;}
																						}#foreach $allele
																						if($A_length == $increment){
																								push @addAllSequences, "1";
																								push @A_file,"\t\tDiagnostic:\t".$key."-->";
																								foreach my $allele(@majorAlleles){push @A_file, $allele." ";}
																								push @A_file, "\n";
																								@A_file	= mhc_associatedAlleles(	Allele => \@majorAlleles,
																																	File_In => \@A_file, 
																																	Sample_Array => \@sampleArray, 
																																	Type => "Major");
																								my @minAlleles	= mhc_get_minorAlleles(	fileIn=> "haploData",
																																			KeyIn=> $key,
																																			HashIn=> \%haploHash);
																								push @A_file, "\t\t\t\t\tMinor Alleles ";
																								foreach my $allele(@minAlleles){push @A_file, $allele." ";}
																								push @A_file, "\n";
																								@A_file	= mhc_associatedAlleles( Allele => \@minAlleles,
																																	File_In	=> \@A_file, 
																																	Sample_Array => \@sampleArray, 
																																	Type => "Minor");
																						}#$A_length==$increment
																				}#A_length
																				else{
																						push @addAllSequences, "1";
																						push @A_file,"\t\tDiagnostic:\t".$key."-->";
																						foreach my $allele(@majorAlleles){push @A_file, $allele." ";}
																						push @A_file, "\n";
																						@A_file = mhc_associatedAlleles(	Allele => \@majorAlleles,
																													File_In => \@A_file, 
																													Sample_Array => \@sampleArray,
																													Type => "Major");
																						my @minAlleles	= mhc_get_minorAlleles(	fileIn=> "haploData",
																															KeyIn=> $key,
																																HashIn=> \%haploHash);
																						push @A_file, "\t\t\t\t\tMinor Alleles ";
																						foreach my $allele4(@minAlleles){push @A_file, $allele4." ";}
																						push @A_file, "\n";
																						@A_file = mhc_associatedAlleles(	Allele	=> \@minAlleles,
																										File_In	=> \@A_file, 
																										Sample_Array => \@sampleArray, 
																										Type => "Minor");
																				}#else
																				my $passed = @addAllSequences;
																				if($passed ==0){}
																}#presentTermArray
														}#arrrayLength
														else{
																foreach my $key(keys %{$haploHash{haploData}}){
																		if($haploHash{haploData}{$key}{'Diagnostic Major'} eq $term){
																				push @A_file, "\t\tDiagnostic:\t".$key."\n";
																				my @minAlleles	= mhc_get_minorAlleles(	fileIn=> "haploData",
																													KeyIn=> $key,
																													HashIn=> \%haploHash);
																				push @A_file, "\t\t\t\t\tMinor Alleles ";
																				foreach my $allele4(@minAlleles){push @A_file, $allele4." ";}
																				push @A_file, "\n";
																				@A_file	= mhc_associatedAlleles(	Allele	=> \@minAlleles,
																											File_In	=> \@A_file, 
																											Sample_Array => \@sampleArray, 
																											Type => "Minor");
																		}#if #ERIKA EDITED LATER
																}#key
														}#else
													#} ERIKA DELTED LATER
												}#$noRepeat
												else{}
								}#searchTermArray
						}#reverseSortedArray
				}#A1Array


#************************ @BArray ********************************

				my @BArray	= grep { $_ =~ /B\_/ } @sampleArray;
				if(@BArray){
						my @revSortedBArray = map{$_->[1]} sort{$b->[0]<=>$a->[0]} map{[substr($_,0,index($_,"\t")),$_]}@BArray;
						push @B_file,"Total Alleles: ".@revSortedBArray."\n";
						push @B2_file,"Total Alleles: ".@revSortedBArray."\n";
						my @usedTerms;

						my $haplotype=0;
						foreach my $element (@revSortedBArray){
								$haplotype++;
								my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $element;
								push @B_file, "H".$haplotype."\t".$allele."\t".$readCount."\n";
								push @B2_file, "H".$haplotype."\t".$allele."\t".$readCount."\n";
								my @haploDataArray	= mhc_allelesToSearch(	Term => "B",
																				Allele => $allele);
								@haploDataArray = uniqueArray(@haploDataArray);
								foreach my $term (@haploDataArray){
										my $noRepeat = 0;
										foreach my $used(@usedTerms){#Has the term previously been used in the current sample
												if	($used eq $term){$noRepeat = 1;}
										}
										if ($noRepeat == 0){
												push @usedTerms, $term;
												my @presentTermArray;
												foreach my $key (keys %{$haploHashB{haploDataB}}){#Get all alleles and suballeles identified 
														if ($haploHashB{haploDataB}{$key}{'Diagnostic Major'} eq $term){
																push @presentTermArray, $key
														}#if
												}#$key
												my $arrayLength = @presentTermArray;#check for multiple alleles
												my @addAllSequences;
												if	($arrayLength > 1){#sub-alleles a,b,c etc
														foreach my $key	(@presentTermArray){
																my $increment= 0;
																my @majorAlleles = mhc_get_majorAlleles(	fileIn=> "haploDataB",
																											KeyIn	=> $key,
																											HashIn=> \%haploHashB);
																my @minAlleles = mhc_get_minorAlleles(	fileIn	=> "haploDataB",
																											KeyIn => $key,
																											HashIn	=> \%haploHashB);
#************************* B2 Allele *************************************************************************
																push @B2_file, "\t\t\tDiagnostic:\t".$key."--> ";
																foreach my $majorAllele	(@majorAlleles){
																		push @B2_file,$majorAllele." ";
																}#foreach my $allele
																push @B2_file, "\n";
																@B2_file = mhc_associatedAlleles(	Allele	=> \@majorAlleles,
																										File_In	=> \@B2_file, 
																										Sample_Array => \@sampleArray,
																										Type => "Major");

																push @B2_file, "\t\t\t\t\tMinor Alleles ";
																foreach my $minAllele	(@minAlleles){push @B2_file, $minAllele." ";}
																push@B2_file, "\n";
																@B2_file = mhc_associatedAlleles(	Allele	=> \@minAlleles,
																										File_In	=> \@B2_file, 
																										Sample_Array => \@sampleArray, 
																										Type => "Minor");
#************************* B Allele *********************************************************************
																my $B_length = @majorAlleles;
																foreach my $majorAllele (@majorAlleles){
																		$majorAllele =~ s/\*/\_/g;
																		my @assMajAllele = grep { $_ =~ /$majorAllele/ } @sampleArray;
																		if	(@assMajAllele){$increment++;}
																}#foreach $allele
																my @tempB_file;
																if($B_length == $increment){#All associated alleles of the sub allele have been found
																		my @haplosIdentified;
																		push @addAllSequences, "1";
																		push @tempB_file, "\tDiagnostic:\t".$key."--> ";
																		foreach my $majorAllele(@majorAlleles){push @tempB_file, $majorAllele." ";}
																		push @tempB_file, "\n";
																		foreach my $majorAllele(@majorAlleles){
																				my @assMajAllele	= grep { $_ =~ /$majorAllele/ } @sampleArray;
																				foreach my $majors(@assMajAllele){
																						for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $majors){push @haplosIdentified, $i+1;}}
																						my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $majors;
																						push @tempB_file, "\t\t\t\t\tMajor:\t".$majorAllele."\t".$readCount."\t".$allele."\n";
																				}
																		}#$allele in @majorAlleles
																		push @tempB_file, "\t\t\t\t\t\tMinor Alleles:\t";
																		foreach my $minAllele(@minAlleles){push @tempB_file, $minAllele." ";}
																		push @tempB_file, "\n";
																		foreach my $minAllele(@minAlleles){
																				$minAllele =~ s/\*/\_/g;
																				my @assMinAllele	= grep { $_ =~ /$minAllele/ } @sampleArray;
																				foreach my $minors(@assMinAllele){
																						for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $minors){push @haplosIdentified, $i+1;}}
																						my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $minors;
																						push @tempB_file, "\t\t\t\t\t\t\tMinor:\t".$minAllele."\t".$readCount."\t".$allele."\n";
																				}
																		}
																		push @B_file, "\t\t\t";
																		for my $haplos(@haplosIdentified){push @B_file, "H".$haplos." ";}
																		push @B_file, @tempB_file;
																}#B_length = $increment
														}# $key in @presentTermArray
												}#if $arrayLength>1
												else{#no sub alleles. Only one Diagnostic Major
														foreach my $key(keys %{$haploHashB{haploDataB}}){
																if($haploHashB{haploDataB}{$key}{'Diagnostic Major'} eq $term){
																		my @tempB_file;
																		my @haplosIdentified;
																		push @addAllSequences, "1";
																		push @tempB_file, "\tDiagnostic:\t".$key."--> ";
																		push @B2_file, "\t\t\tDiagnostic:\t".$key."--> ";
																			my @majorAlleles = mhc_get_majorAlleles(	fileIn=> "haploDataB",
																														KeyIn => $key,
																														HashIn=> \%haploHashB);
																		my @minAlleles = mhc_get_minorAlleles(	fileIn => "haploDataB",
																													Key => $key,
																													HashIn => \%haploHashB);
																		foreach my $majorAllele(@majorAlleles){push @tempB_file, $majorAllele." ";push @B2_file, $majorAllele." ";}
																		push @tempB_file, "\n";
																		push @B2_file, "\n";
																		foreach my $majorAllele(@majorAlleles){
																				$majorAllele=~ s/\*/\_/g;
																				my @assMajAllele	= grep { $_ =~ /$majorAllele/ } @sampleArray;
																				foreach my $majors(@assMajAllele){
																						for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $majors){push @haplosIdentified, $i+1;}}
																						my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $majors;
																						push @tempB_file, "\t\t\t\t\tMajor:\t".$majorAllele."\t".$readCount."\t".$allele."\n";
																				}
																		}
																		@B2_file = mhc_associatedAlleles(	Allele => \@majorAlleles,
																												File_In => \@B2_file, 
																												Sample_Array => \@sampleArray, 
																												Type => "Major");

																		push @tempB_file, "\t\t\t\t\t\tMinor Alleles\t ";
																		push @B2_file, "\t\t\t\t\tMinor Alleles ";
																		foreach my $minAllele(@minAlleles){push @tempB_file, $minAllele." ";push @B2_file, $minAllele." ";}
																		push @tempB_file, "\n";
																		push @B2_file, "\n";
																		foreach my $minAllele(@minAlleles){
																				$minAllele =~ s/\*/\_/g;
																				my @assMinAllele = grep { $_ =~ /$minAllele/ } @sampleArray;
																				foreach my $minors(@assMinAllele){
																						for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $minors){push @haplosIdentified, $i+1;}}
																						my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $minors;
																						push @tempB_file, "\t\t\t\t\t\t\tMinor:\t".$minAllele."\t".$readCount."\t".$allele."\n";
																				}
																		}
																		@B2_file = mhc_associatedAlleles(	Allele => \@minAlleles,
																										File_In => \@B2_file, 
																										Sample_Array => \@sampleArray, 
																										Type => "Minor");
																		push @B_file, "\t\t\t";
																		for my $haplos(@haplosIdentified){push @B_file, "H".$haplos." ";}
																		push @B_file, @tempB_file;
																}#if Diagnostic Major
														}#key
												}#else: $arrayLength=<1
												my $passed = @addAllSequences;
												if($passed ==0){
														foreach my $key(keys %{$haploHashB{haploDataB}}){
																if($haploHashB{haploDataB}{$key}{'Diagnostic Major'} eq $term){
																		my @tempB_file;
																		my @haplosIdentified;
																		push @tempB_file, "\tDiagnostic Incomplete:\t".$key."--> ";
																		my @majorAlleles = mhc_get_majorAlleles(	fileIn => "haploDataB",
																													KeyIn => $key,
																													HashIn=> \%haploHashB);
																		my @minAlleles = mhc_get_minorAlleles( fileIn	=> "haploDataB",
																												  KeyIn		=> $key,
																												  HashIn	=> \%haploHashB);
																		foreach my $majorAllele(@majorAlleles){push @tempB_file, $majorAllele." ";}
																		push @tempB_file, "\n";
																		foreach my $majorAllele(@majorAlleles){
																		$majorAllele =~ s/\*/\_/g;
																		# print $majorAllele."\n";
																		my @assMajAllele	= grep { $_ =~ /$majorAllele/ } @sampleArray;
																		foreach my $majors(@assMajAllele){
																				for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $majors){push @haplosIdentified, $i+1;}}
																				my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $majors;
																				push @tempB_file, "\t\t\t\t\tMajor:\t".$majorAllele."\t".$readCount."\t".$allele."\n";}
																		}
																		push @tempB_file, "\t\t\t\t\t\t\tMinor Alleles ";
																		foreach my $minAllele	(@minAlleles){push @tempB_file, $minAllele." ";}
																		push @tempB_file, "\n";
																		foreach my $minAllele(@minAlleles){
																				$minAllele 			=~ s/\*/\_/g;
																				my @assMinAllele	= grep { $_ =~ /$minAllele/ } @sampleArray;
																				foreach my $minors(@assMinAllele){
																				for (my $i = 0;$i<@revSortedBArray;$i++){if($revSortedBArray[$i] eq $minors){push @haplosIdentified, $i+1;}}
																						my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $minors;
																						push @tempB_file, "\t\t\t\t\t\t\t\tMinor:\t".$minAllele."\t".$readCount."\t".$allele."\n";}
																		}#foreach
																		push @B_file, "\t\t\t";
																		for my $haplos(@haplosIdentified){push @B_file, "H".$haplos." ";}
																		push @B_file, @tempB_file;
																}#if DiagMaj = term
														}#foreach haploDataB
												}#if passed
										}#noRepeat
								}#term
						}#reverseSortedArray
				}#@BArray
		}#Each Sample



#****************************** Data ****************************************
push @data_file, "Allele\t";
foreach my $sample(@uniqueSampleArray){push @data_file, $sample."\t";}
push @data_file, "Representation\t%\n";
		foreach my $uniqueAllele(@uniqueAlleleArray){
				my $alleleInSamples 			=0;
				push @data_file, $uniqueAllele."\t";
				my $alleleFound	=0;
				for	(my $i=0; $i<@uniqueSampleArray; $i++){
						my @sampleArray;
						foreach my $line (@fileArray){
								my($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $line;
								if 	($sample eq $uniqueSampleArray[$i]){ push @sampleArray, $line; }
						}
						my @statsArray = @sampleArray;
						my @revStats = map{$_->[1]} sort{$b->[0]<=>$a->[0]} map{[substr($_,0,index($_,"\t")),$_]}@statsArray;
						for(my $j = 0; $j <@revStats; $j++){
								my	($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $revStats[$j];
								my $statsAllele		= $allele;
								if($statsAllele eq $uniqueAllele){
										$alleleInSamples++;
										$alleleFound 	=1;
										my $rank = $j+1;
										push @data_file, "(".$rank.")".$readCount;
										$j = @revStats;
								}
						}
						push @data_file, "\t";
				}
				push @data_file, $alleleInSamples."/".$uniqueLength."\t".($alleleInSamples/$uniqueLength);
				push @data_file, "\n";
		}

		write_file			($var_pat."/_reports/rawData.txt",@fileArray);
		write_file			($var_pat."/_reports/data.txt",@data_file);
		write_file			($var_pat."/_reports/A.txt",@A_file);
		write_file			($var_pat."/_reports/B.txt",@B_file);
		write_file			($var_pat."/_reports/B_detail.txt",@B2_file);

		my (@align_excel)			= <$var_pat/_reports/*.txt>; 
		foreach(@align_excel){ 
				$_= "_reports/".returnFile($_);
		}
		ExcelfromArray(	Worksheets => \@align_excel,
									Path => $var_pat,
									Name=> "projectSummary.xlsx"); 
}#mhc_A_B

sub mhc_allelesToSearch{
		my (%args)			=(@_);
		my $searchTerm = $args{Term};
		my $allele = $args{Allele};
		my @haploDataArray;
		my $haploDataSearch;
		$allele =~ s/\|/\&\_/g;
		my @sampleSplit		 			= split /\_/, $allele;
		for (my $index=0;$index<@sampleSplit;$index++){
				if ($sampleSplit[$index] eq $searchTerm){
						my $alleleNumber= $sampleSplit[$index+1];
						my $fragment;
						if ($alleleNumber =~ /g/){$fragment = substr $alleleNumber ,0,index($alleleNumber,"g"); }
						else{$fragment = $alleleNumber;}
						$haploDataSearch= $searchTerm."*".$fragment;
						push  @haploDataArray, $haploDataSearch;
				}#if
		}#for
		@haploDataArray	= uniqueArray(@haploDataArray);
		return @haploDataArray;
}

sub mhc_get_minorAlleles{
		my (%args) = (@_);
		my $fileNam = $args{fileIn};
		my $key	= $args{KeyIn};
		my %haploHash		= %{$args{HashIn}};
		my @minAlleles;
		if($haploHash{$fileName}{$key}{'Minor 1'}){push minAlleles, $haploHash{$fileName}{$key}{'Minor 1'};}
		if($haploHash{$fileName}{$key}{'Minor 2'}){push minAlleles, $haploHash{$fileName}{$key}{'Minor 2'};}
		if($haploHash{$fileName}{$key}{'Minor 3'}){push minAlleles, $haploHash{$fileName}{$key}{'Minor 3'};}
		if($haploHash{$fileName}{$key}{'Minor 4'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 4'};}
		if($haploHash{$fileName}{$key}{'Minor 5'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 5'};}
		if($haploHash{$fileName}{$key}{'Minor 6'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 6'};}
		if($haploHash{$fileName}{$key}{'Minor 7'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 7'};}
		if($haploHash{$fileName}{$key}{'Minor 8'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 8'};}
		if($haploHash{$fileName}{$key}{'Minor 9'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 9'};}
		if($haploHash{$fileName}{$key}{'Minor 10'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 10'};}
		if($haploHash{$fileName}{$key}{'Minor 11'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 11'};}
		if($haploHash{$fileName}{$key}{'Minor 12'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 12'};}
		if($haploHash{$fileName}{$key}{'Minor 13'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 13'};}
		if($haploHash{$fileName}{$key}{'Minor 14'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 14'};}
		if($haploHash{$fileName}{$key}{'Minor 15'}){push @minAlleles, $haploHash{$fileName}{$key}{'Minor 15'};}

@minAlleles 		= grep {$_ !~ /\-/} @minAlleles;

return	@minAlleles;
}

sub mhc_get_majorAlleles{
		my (%args) = (@_);
		my $fileName = $args{fileIn};
		my $key	= $args{KeyIn};
		my %haploHash = %{$args{HashIn}};
		my @majAlleles;

		if($haploHash{$fileName}{$key}{'Major 2'}){push 		@majAlleles, $haploHash{$fileName}{$key}{'Major 2'};}
		if($haploHash{$fileName}{$key}{'Major 3'}){push 		@majAlleles, $haploHash{$fileName}{$key}{'Major 3'};}
		if($haploHash{$fileName}{$key}{'Major 4'}){push 		@majAlleles, $haploHash{$fileName}{$key}{'Major 4'};}
		if($haploHash{$fileName}{$key}{'Major 5'}){push 		@majAlleles, $haploHash{$fileName}{$key}{'Major 5'};}
		if($haploHash{$fileName}{$key}{'Major 6'}){push 		@majAlleles, $haploHash{$fileName}{$key}{'Major 6'};}
		@majAlleles = grep {$_ !~ /\-/} @majAlleles;
		
		return @majAlleles;
}

sub mhc_associatedAlleles{
		my (%args)	=(@_);
		my @majorAlleles = (@{$args{Allele}});
		my @file=(@{$args{File_In}});
		if (@majorAlleles){
				foreach my $value (@majorAlleles){
						$value =~ s/\*/\_/g;
						my @assocAllele= grep { $_ =~ /$value/ } (@{$args{Sample_Array}});
						foreach my $associated (@assocAllele){
							my ($readCount, $sample, $allele, $percentReads, $mappedReads, $totalReads) = split "\t", $associated;
							if($args{Type} eq "Minor"){push 			@file,"\t\t";}
							push @file, "\t\t\t\t".$args{Type}.":\t".$value."\t".$readCount."\t".$allele."\n";
						}#foreach
				}#foreach
		}#if
		return @file;
}
