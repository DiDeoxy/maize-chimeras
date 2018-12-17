#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::SeqFeature::Store;

my $usage = "Usage: $0 snps.txt out_file.txt";

my $snps = shift @ARGV || die $usage;
my $out_file = shift @ARGV || die $usage;

my $gff_db = Bio::DB::SeqFeature::Store->new(	-adpator => 'DBI::mysql',
												-dsn => "DBI:mysql:database=gff;host=localhost;user=max",
												-index_subfeatures => 0);

my %id_postns;
open (IN, $snps);
while (my $snp = <IN>) {
	my ($id, $postn, $base) = split('\s', $snp);
	push (@{ $id_postns{$id} }, [$postn, $base]);
}
close IN;

my @types = qw(exon);

open (OUT, ">$out_file") || die "Failed to open $out_file";
for my $id (keys %id_postns) {
	my @features = $gff_db->get_features_by_name($id);
	
	for my $feature (@features) {
		my @subfeatures = $feature->get_SeqFeatures(@types);
		
		for my $postn_base (@{ $id_postns{$id} }) {
			my ($postn, $base) = @$postn_base;
			my $postn_dif = $postn;
			
			if ($feature->strand == 1) {
				
				for my $subfeature (@subfeatures) {
					my $exact_postn = $subfeature->start + $postn_dif - 1;
					if ($exact_postn <= $subfeature->end) {
						print OUT $feature->seq_id, "\t", $feature->display_name, "\t$exact_postn\t$base\t", $feature->strand, "\n";
						last;
					} else {
						$postn_dif = ($subfeature->start + $postn_dif - $subfeature->end) - 1; 
					}
				}
				
			} else {
				
				for my $subfeature (@subfeatures) {
					my $exact_postn = ($subfeature->end - $postn_dif) + 1;
					if ($exact_postn >= $subfeature->start) {
						print OUT $feature->seq_id, "\t", $feature->display_name, "\t$exact_postn\t$base\t", $feature->strand, "\n";
						last;
					} else {
						$postn_dif = ($subfeature->start - ($subfeature->end - $postn_dif)) - 1; 
					}
				}
				
			}
		}
	}
}
close OUT;