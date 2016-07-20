use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];
my $start_run = time();

open (IN, "<$inputfile");
open (OUT, ">$outputfile");

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous',
	-port => 3337,
	);

# Connect to the databases:
my $variation_adaptor = $registry->get_adaptor('homo_sapiens', 'variation', 'variation');
my $ldfc_adaptor = $registry->get_adaptor('homo_sapiens', 'variation', 'ldfeaturecontainer');
$ldfc_adaptor->max_snp_distance(500_000);
my $pop_adaptor = $registry->get_adaptor('homo_sapiens', 'variation', 'population');
$variation_adaptor->db->use_vcf(1);

my $ld_population = $pop_adaptor->fetch_by_name('1000GENOMES:phase_3:CEU');

print OUT "chrnum\tpos_ldSNP\tldSNP\tpos_tagSNP\ttagSNP\tr2\n";

# Loop through all SNPs available and find SNPs in LD
while(<IN>) {
	chomp;
	my $variation_name = $_;
	my $variation = $variation_adaptor->fetch_by_name($variation_name);
	my @var_features;

	if ($variation) {
		@var_features = @{ $variation->get_all_VariationFeatures() };
	} else {
		print 'failing variation name: ', $variation_name, "\n";
		next;
	}

	foreach my $vf (@var_features) {
		if ($vf->slice->is_reference) {

			my $rsid = $vf->name;
			print $rsid, "\n";
			my $start = $vf->start;
			my $region = $vf->seq_region_name;
			
			my $ldfc = $ldfc_adaptor->fetch_by_VariationFeature($vf, $ld_population);
			my @ld_values = @{ $ldfc->get_all_ld_values() };
	      	foreach my $ld_hash (@ld_values) {
	        	my $r2 = $ld_hash->{r2};
		    	my $variation_name1 = $ld_hash->{variation_name1};
		    	my $variation_name2 = $ld_hash->{variation_name2};
		    	my $pos1 = $ld_hash->{variation1}->seq_region_name();
		    	my $pos2 = $ld_hash->{variation2}->seq_region_name();
		    	my $start1 = $ld_hash->{variation1}->start;
		    	my $start2 = $ld_hash->{variation2}->start;

		    	if ($variation_name1 eq $rsid) {
		    		print OUT "$pos2\t$start2\t$variation_name2\t$start1\t$variation_name1\t$r2\n";
		    	} else {
		    		print OUT "$pos1\t$start1\t$variation_name1\t$start2\t$variation_name2\t$r2\n";
		    	}
		    }
      	} else {
      		print "Failed to get, not a reference allele: ", $vf->name, "\n";
      	}
	}
}

close OUT;
close IN;
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";