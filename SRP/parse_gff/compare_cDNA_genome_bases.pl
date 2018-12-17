#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Fasta;

my $usage = "Usage: $0 genome.fasta genome_postns_info.txt disagreeing_bases.txt tophat_out.txt";

my $genome_fasta = shift @ARGV || die $usage;
my $genome_postns = shift @ARGV || die $usage;
my $disagree_out = shift @ARGV || die $usage;
my $tophat_out = shift @ARGV || die $usage;

my $db = Bio::DB::Fasta->new("$genome_fasta");
my @ids = $db->get_all_primary_ids;

my %chromosomes;
for my $chromosome (@ids) {
	$chromosome =~ /chromosome:AGPv2:(\w+)/;
	
	if ($1 eq "mitochondrion") {
		$chromosomes{"Mt"} = $chromosome;
	} else {
		$chromosomes{$1} = $chromosome;
	}
}

my %chr_cdna_strand_info;
open (IN, $genome_postns) || die "Failed at open $genome_postns";
while (my $line = <IN>) {
	my ($chr, $cdna, $genome_postn, $cdna_base, $strand) = split '\s', $line;
	
	push @{ $chr_cdna_strand_info{"$chr,$cdna,$strand"} }, [$genome_postn, $cdna_base];
}
close IN;

open OUT1, ">$disagree_out" || die "Failed at open $disagree_out";
open OUT2, ">$tophat_out" || die "Failed at open $tophat_out";
foreach my $chr_cdna_strand (keys %chr_cdna_strand_info) {
	
	my ($chr, $cdna, $strand) = split ',', $chr_cdna_strand;
	
	my $chromosome = $chromosomes{$chr};
	
	for my $info (@{ $chr_cdna_strand_info{$chr_cdna_strand} }) {
		
		my ($genome_postn, $cdna_base) = @$info;
		
		my $genome_base = $db->seq("$chromosome", $genome_postn => $genome_postn);
		
		if ($strand == 1) {
			if ($cdna_base ne $genome_base) {
				print OUT1 "$cdna_base\t$genome_base\n";
			}
		} else {
			$genome_base =~ tr/ACTG/TGAC/;
			if ($cdna_base ne $genome_base) {
				print OUT1 "$cdna_base\t$genome_base\n";
			}
		}
		
		print OUT2 "$chromosome\t$genome_postn\n";
		
	}
}
close OUT1;
close OUT2;