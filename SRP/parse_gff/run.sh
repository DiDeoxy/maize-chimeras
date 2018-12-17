# ## Builds a MySQL database from the gff file
# perl build_db_from_gff.pl ~/Downloads/ZmB73_5b_FGS.gff.gz
#
# ## Maps ZmB73 FGS cDNA mismatch positions to the full genome
# perl map_snp_positions_to_genome.pl ~/Downloads/Z_b+f_snp_postns.txt ~/Downloads/genomic_postns_info.txt
#
# ## Checks that the genome positions and the cDNA positions return the same base
perl compare_cDNA_genome_bases.pl ~/Downloads/genome ~/Downloads/genomic_postns_info.txt ~/Downloads/disagreeing_bases.txt ~/Downloads/genomic_postns.txt

						


