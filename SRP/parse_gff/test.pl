#!/usr/bin/perl
use strict;
use warnings;

my @relative_pos = 1..100;
my @relative_pos = reverse @relative_pos;
my @actual_pos = (100..109,125..174,210..249);

for (my $i = 0; $i < @relative_pos; $i++) {
	print "$relative_pos[$i]\t$actual_pos[$i]\n";
}