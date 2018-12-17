#!/usr/bin/perl
use strict;
use warnings;
use Bio::SearchIO;

#Max H., June 1, adapted from blast_sear_io2.pl by Ann Meyer
# this pareses full blast outputs and prints to file useful summary statistics

my $usage = "perl script.pl X_on_Ydb.txt outfile.txt";

my $blast_file = shift || die $usage;
my $out_file = shift || die $usage;

#now read in the blast file
my $searchIO = new Bio::SearchIO(-format => 'blast',
						   -file   => $blast_file);
			 
open OUT, ">$out_file" || die $usage;

while (my $result = $searchIO->next_result) {
	while (my $hit = $result->next_hit) {
		while (my $hsp = $hit->next_hsp) {
			#only want hsp with alignmnet lengths greater than 100 and percent identity > 75
			if ($hsp->length('total') >= 100 && $hsp->percent_identity >= 75) {
			
				#put all the mismatches in query strand into an array
				my @quer_mis;
				foreach ($hsp->seq_inds('query','nomatch') ) {
					push (@quer_mis, $_);
				}
				my $quer_mis_string = join(',', @quer_mis);
				unless (length $quer_mis_string) {
					$quer_mis_string = "NA";
				}

				#put all the mismatches in subject strand into an array
				my @sub_mis;
				foreach ($hsp->seq_inds('hit','nomatch') ) {
					push (@sub_mis, $_);
				}
				my $sub_mis_string = join(',', @sub_mis);
				unless (length $sub_mis_string) {
					$sub_mis_string = "NA";
				}

				#put query gaps in array
				my @gap_quer_array = $hsp->seq_inds('query'=>'gaps');
				my $gap_pos_quer = join(',', @gap_quer_array);
				unless (ength $gap_pos_quer) {
					$gap_pos_quer = "NA";
				}
				
				#put query gaps in array
				my @gap_sub_array = $hsp->seq_inds('hit'=>'gaps');
				my $gap_pos_sub = join(',', @gap_sub_array);
				unless (length $gap_pos_sub) {
					$gap_pos_sub = "NA";
				}
				
				#get query range
				my $quer_range = join(',', $hsp->range('query'));
				
				#get sub range
				my $sub_range = join(',', $hsp->range('hit'));


				print OUT $result->query_name,"\t",$hit->name,"\t",$quer_range,"\t",$sub_range,"\t",$hsp->strand('query'),"\t";
				print OUT $hsp->strand('hit'),"\t",$hsp->percent_identity,"\t",$hsp->evalue,"\t",$hsp->length('total'),"\t";
				print OUT $hsp->gaps,"\t",$quer_mis_string,"\t",$sub_mis_string,"\t",$gap_pos_quer,"\t",$gap_pos_sub,"\n";
			}
		}
	}	
}
close OUT;