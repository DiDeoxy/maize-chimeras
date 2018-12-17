#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

my $usage = "Usage: $0 file.gff";

my $gff = shift @ARGV || die $usage;

my $gff_db = Bio::DB::SeqFeature::Store->new(	-adpator => 'DBI::mysql',
												-dsn => "DBI:mysql:database=gff;host=localhost;user=max",
												-index_subfeatures => 0,
												-create => 1);

my $gff_loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(	-store => $gff_db,
																-verbose => 1,
																-fast => 1);

$gff_loader->load($gff);