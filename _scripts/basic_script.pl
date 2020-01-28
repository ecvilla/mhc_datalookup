#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use Data::Dumper;



sub directoryStructure {
my(%args) = (@_);
foreach my $folder(@{$args{Folders}}){
unless (-d $args{Path}."/".$folder) { mkdir $args{Path}."/".$folder; }
}
}

sub uniqueArray {
my %h;
return grep { !$h{$_}[0]++ } @_;
}

sub returnFile {
my ($path)= (@_);
my (@elements) = split "/", $path;
return $elements[@elements-1];
}

sub ExcelfromArray {
		my(%args) = (@_);
		my ($workbook)= Excel::Writer::XLSX -> new($args{Path}."/".$args{Name});
		my $format= $workbook -> add_format();
		$format-> set_align('left'); 
		$format-> set_font('Courier New'); 
		my (@sheet_array)= @{$args{Worksheets}};

		foreach my $sheet(@sheet_array){
				my @lines= read_file($args{Path}."/".$sheet);
				my $worksheet= $workbook -> add_worksheet(bk_basic_removeFileType(bk_basic_returnFile($sheet)));
				for(my $row=0; $row<@lines; $row++){
						my @values= split /\t/, $lines[$row];
						for(my $col=0; $col<@values; $col++){
								$worksheet-> write($row, $col, $values[$col], $format);
						}
				}
		}
}