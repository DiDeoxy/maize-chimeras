#!/usr/bin/perl
use strict;
use warnings;

# Max H., June 2, 2015.
# This script finds SNP and indel positions in default format blast output parsed by parse_blast_searchIO.pl
# Recommended but not required is filtering for RBHs using get_rbh_v2.pl 
# and finding putative alleles by aligning to an anchoring genome (get_putative_alleles_BM_BZ_MZ_v2.pl)

my $usage = "perl script.pl putative_alleles.txt out_snp_pos.txt out_query_gap_indels.txt out_subject_gap_indels.txt \n";

my $putative_alleles = shift || die $usage;
my $out_file1 = shift || die $usage;
my $out_file2 = shift || die $usage;
my $out_file3 = shift || die $usage;

my %blast_hits = &parser($putative_alleles);

my $query_total = 0;
my $subject_total = 0;
# for printing out SNP positions
open OUT1, ">$out_file1" || die $usage;
# for printing out indels where the query sequence has a gap
open OUT2, ">$out_file2" || die $usage;
# for printing out indels where the subject sequence has a gap
open OUT3, ">$out_file3" || die $usage;
# looks at each query-subject pair
for my $query_id (keys %blast_hits) {
	for my $subject_id (keys %{ $blast_hits{$query_id} }) {
		# making the query-subject pairs info available
		my ($query_range, $subject_range, $query_strand, $subject_strand, $num_gaps, $query_mismatch_postns, 
		$subject_mismatch_postns, $query_gap_postns, $subject_gap_postns) = split('\s', $blast_hits{$query_id}{$subject_id});

		# arrays for holding mismatch and gap positions
		my @Query_mismatch_postns;
		my @Query_gap_postns;
		# array for holding both mismatch and gap positions
		my @Query_postns;
		my $num_query_mismatch = 0;
		my $num_query_gap = 0;
		
		# when there are query mismatches fills the appropriate array
		if ($query_mismatch_postns ne "NA") {
			@Query_mismatch_postns = split(',', $query_mismatch_postns);
			$num_query_mismatch = ($#Query_mismatch_postns + 1);
		}
		
		# when query gaps are present fills the appropriate array
		if ($query_gap_postns ne "NA") {
			# sorts the array so that its elements values are in ascending order
			# this is important as it allows mismatches or indels of the alignment to be identified
			# by moving from pairs of such features between the subject and query strand in a consecutive manner
			@Query_gap_postns = sort {$a <=> $b} split(',', $query_gap_postns);
			$num_query_gap = ($#Query_gap_postns + 1);
		}
		
		# fills the query positions array in a manner dependent upon whether there are gaps, mismatches, or both
		# this is done so that we can eliminate indels by comparing this new array to the gap array
		# mismatch positions will not have  identical positions in the gap array (or at least one less)
		if ($query_mismatch_postns ne "NA" && $query_gap_postns ne "NA") {
			@Query_postns = sort {$a <=> $b} (@Query_mismatch_postns, @Query_gap_postns);
		} elsif ($query_mismatch_postns ne "NA") {
			@Query_postns = sort {$a <=> $b} @Query_mismatch_postns;
		} elsif ($query_gap_postns ne "NA") {
			@Query_postns = sort {$a <=> $b} @Query_gap_postns;
		}
		
		# same as above but for the subject array
		my @Subject_mismatch_postns;
		my @Subject_gap_postns;
		my @Subject_postns;
		
		my $num_subject_mismatch = 0;
		my $num_subject_gap = 0;
		
		# the subject sequence can be aligned in the forward or backwards direction, 1 indicates forward
		# in this case the mismatches and gaps are parsed identically to the query strand (which is always forward)
		if ($subject_strand == 1) {
			
			if ($subject_mismatch_postns ne "NA") {
				@Subject_mismatch_postns = split(',', $subject_mismatch_postns);
				$num_subject_mismatch = ($#Subject_mismatch_postns + 1);
			}
			
			if ($subject_gap_postns ne "NA") {
				@Subject_gap_postns = sort {$a <=> $b} split(',', $subject_gap_postns);
				$num_subject_gap = ($#Subject_gap_postns + 1);
			}
			
			if ($subject_mismatch_postns ne "NA" && $subject_gap_postns ne "NA") {
				@Subject_postns = sort {$a <=> $b} (@Subject_mismatch_postns, @Subject_gap_postns);
			} elsif ($subject_mismatch_postns ne "NA") {
				@Subject_postns = sort {$a <=> $b} @Subject_mismatch_postns;
			} elsif ($subject_gap_postns ne "NA") {
				@Subject_postns = sort {$a <=> $b} @Subject_gap_postns;
			}
		
		# this is for when the subject strand is aligned backwards
		# in this case the subject alignment features must be ordered in descending order
		} else {
			
			if ($subject_mismatch_postns ne "NA") {
				@Subject_mismatch_postns = split(',', $subject_mismatch_postns);
				$num_subject_mismatch = ($#Subject_mismatch_postns + 1);
			}
			
			if ($subject_gap_postns ne "NA") {
				@Subject_gap_postns = sort {$b <=> $a} split(',', $subject_gap_postns);
				$num_subject_gap = ($#Subject_gap_postns + 1);
			}
			
			if ($subject_mismatch_postns ne "NA" && $subject_gap_postns ne "NA") {
				@Subject_postns = sort {$b <=> $a} (@Subject_mismatch_postns, @Subject_gap_postns);
			} elsif ($subject_mismatch_postns ne "NA") {
				@Subject_postns = sort {$b <=> $a} @Subject_mismatch_postns;
			} elsif ($subject_gap_postns ne "NA") {
				@Subject_postns = sort {$b <=> $a} @Subject_gap_postns;
			}
			
		}
		
		# counts the total number of mismatches that the logic should output
		$query_total += $num_query_mismatch - $num_subject_gap;
		$subject_total += $num_subject_mismatch - $num_query_gap;
		
		# index values for which gap value to compare for the appropriate sequence
		my $query_gap_index = 0;
		my $subject_gap_index = 0;
		# both _postn arrays will have the same number of elements therefore we can loop through them as pairs with a for
		for (my $index = 0; $index <= $#Query_postns; $index++) {
			
			# counts the number of concurrent occurrences of the current position value in the proceeding position values inclusive
			my $query_count1 = 0;
			while (	exists $Query_postns[$index + $query_count1] &&
						$Query_postns[$index] == $Query_postns[$index + $query_count1]) {
				$query_count1++;
			}
			
			# same for gap positions
			my $query_count2 = 0;
			while (	exists $Query_gap_postns[$query_gap_index + $query_count2] &&
						$Query_gap_postns[$query_gap_index] == $Query_gap_postns[$query_gap_index + $query_count2]) {
				$query_count2++;
			}
			
			# same for subject
			my $subject_count1 = 0;
			while (	exists $Subject_postns[$index + $subject_count1] &&
						$Subject_postns[$index] == $Subject_postns[$index + $subject_count1]) {
				$subject_count1++;
			}
			
			my $subject_count2 = 0;
			while (	exists $Subject_gap_postns[$subject_gap_index + $subject_count2] &&
						$Subject_gap_postns[$subject_gap_index] == $Subject_gap_postns[$subject_gap_index + $subject_count2]) {
				$subject_count2++;
			}
			
			# uses these counts to find situations where there are more occurrences of the current value in the _postns array than the gap array
			# in this situation the query position is not a gap, but would be printed out as such if treated normally because the current query gap position
			# would be equal to the current _postn position. A second possibility is that there is a subject gap in this pair, therefore we test for this and if not
			# present this is a mismatch and is printed as such. All this only applies when there are more than 0 gaps remaining in the query alignment.
			if ($query_count1 > $query_count2 && $query_count2 != 0 ) {
				#asks if the current subject postn is equivalent to the current gap position, if so we have subject gap indel
				if (exists $Subject_gap_postns[$subject_gap_index] &&
					$Subject_postns[$index] == $Subject_gap_postns[$subject_gap_index]) {
					# prints out subject gap indel	
					print OUT3 "$query_id\t$subject_id\t" . $Query_postns[$index] . "\t$Subject_gap_postns[$subject_gap_index]," . ($Subject_gap_postns[$subject_gap_index] + 1) . "\n";
					# increments the subject index so the next comparison is with the proceeding gap position
					$subject_gap_index++;
					
				} else {
					# when the other options fail it must be a mismatch
					print OUT1 "$query_id\t$subject_id\t" . $Query_postns[$index] . "\t" . $Subject_postns[$index] . "\t" . $query_range . "\t" . $subject_range . "\n";
					
				}
			# same as above but for the subject, therefore everything is inverse	
			} elsif ($subject_count1 > $subject_count2 && $subject_count2 != 0) {
			
				if (exists $Query_gap_postns[$query_gap_index] &&
					$Query_postns[$index] == $Query_gap_postns[$query_gap_index]) {
					
					print OUT2 "$query_id\t$subject_id\t$Query_gap_postns[$query_gap_index]," . ($Query_gap_postns[$query_gap_index] + 1) . "\t" . $Subject_postns[$index] . "\n";
					$query_gap_index++;
					
				} else {
				
					print OUT1 "$query_id\t$subject_id\t" . $Query_postns[$index] . "\t" . $Subject_postns[$index] . "\t" . $query_range . "\t" . $subject_range . "\n";
					
				}
			# all other situations are treated normally as there is no confusion in the remainder as to whether a position is an indel or a mismatch
			} else {
			
				if (exists $Query_gap_postns[$query_gap_index] &&
					$Query_postns[$index] == $Query_gap_postns[$query_gap_index]) {

					print OUT2 "$query_id\t$subject_id\t$Query_gap_postns[$query_gap_index]," . ($Query_gap_postns[$query_gap_index] + 1) . "\t" . $Subject_postns[$index] . "\n";
					$query_gap_index++;
					
				} elsif ( exists $Subject_gap_postns[$subject_gap_index] &&
							$Subject_postns[$index] == $Subject_gap_postns[$subject_gap_index]) {
					
					print OUT3 "$query_id\t$subject_id\t" . $Query_postns[$index] . "\t$Subject_gap_postns[$subject_gap_index]," . ($Subject_gap_postns[$subject_gap_index] + 1) . "\n";
					$subject_gap_index++;

				} else {
				
					print OUT1 "$query_id\t$subject_id\t" . $Query_postns[$index] . "\t" . $Subject_postns[$index] . "\t" . $query_range . "\t" . $subject_range . "\n";
				
				}
			
			}
		}
	}
}
close IN;
close OUT1;
close OUT2;
close OUT3;

# prints out the number of mismatches expected
print "$query_total\t$subject_total\n";

#creates hashes of the each parsed blast line's info
sub parser {
	my $rbhs = shift @_;
	# a hash for holding hashes of top matches
	my %parsed;
	open (IN, $rbhs) || die "Could not open $rbhs for reading.\n";
	while (my $line = <IN>) {
		$line =~ s/\cM//;
		chomp $line;
		
		# taking input line and assigning its components to appropriate scalars
		my ($query_id, $subject_id, $query_range, $subject_range, $query_strand, $subject_strand,
		$percent_identity, $evalue, $align_length, $num_gaps, $query_mismatch_pos, 
		$subject_mismatch_pos, $query_gap_postns, $subject_gap_postns) = split('\s', $line);
		
		$parsed{$query_id}{$subject_id} = "$query_range\t$subject_range\t$query_strand\t$subject_strand\t$num_gaps\t$query_mismatch_pos\t$subject_mismatch_pos\t$query_gap_postns\t$subject_gap_postns";
	}
	close IN;
	return %parsed;
}