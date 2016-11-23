# Calculate mutual exclusive interactors
# Follows the same approach as TCGA mutual exclusivity analyses
# Swati Kaushik March 21 2016

use strict;
use warnings;
use Statistics::R;
use List::MoreUtils qw(uniq);

#read the input file
my $filename = $ARGV[0];
my $mutation_file = $ARGV[1];

open IN, "$filename"  or die $!;
my @interactor_file = <IN>;

my $gene_name = 'XX';

open AH, "$mutation_file" or die $!;
my @tcga_mutation_file =<AH>;

open OUTALL, ">Odds-ratio-ints.txt" or die $!;

foreach my $interactor_name (@interactor_file){
	
	print OUTALL "\n";	
	$interactor_name =~s/\n|\t|\s//g;
	print OUTALL "$interactor_name\t";
	print "$interactor_name\t";
	
	my %XX_mutated_samples; my %interactor_mutated_sample; my %samples_not_altered_G1_G2; my %gene_sample_pairs; my %unique_samples;
	

	#open OUTE, ">XX_mutated_samples" or die $!;
	#open OUTG, ">G2_mutated_samples" or die $!;
	#open OUTB, ">samples_NOT_G1_G2_mutated_samples" or die $!;
	#open OUTC, ">common_mutated_samples" or die $!;

	foreach (@tcga_mutation_file){

		my @mutation_rows = split("\t", $_);
	
		if (($mutation_rows[0] eq $gene_name) && ($mutation_rows[8] ne "Nonsense_Mutation") && ($mutation_rows[8] ne "Silent"))
			{	
				$mutation_rows[15] = substr($mutation_rows[15], 0,15);	
				#print OUTE "$mutation_rows[0]\t$mutation_rows[15]\n";
				$XX_mutated_samples{$mutation_rows[15]}=1;
			}	

		if (($mutation_rows[0] eq $interactor_name) && ($mutation_rows[8] ne "Nonsense_Mutation") && ($mutation_rows[8] ne "Silent"))
			{	
				$mutation_rows[15] = substr($mutation_rows[15], 0,15);	
				#print OUTG "$mutation_rows[0]\t$mutation_rows[15]\n";
				$interactor_mutated_sample{$mutation_rows[15]}=1;
			}	
		
		if (($mutation_rows[8] ne "Nonsense_Mutation") && ($mutation_rows[8] ne "Silent"))
			{	
				$mutation_rows[15] = substr($mutation_rows[15], 0,15);	
				$gene_sample_pairs{$mutation_rows[0]}{$mutation_rows[15]} =0;
				$unique_samples{$mutation_rows[15]} =0;
			}

	}

	#find samples which have both XX and interactor mutations. Compare XX and interactor mutation hash
	my @common_samples=();

	foreach (keys %XX_mutated_samples){

		push (@common_samples, $_) if exists $interactor_mutated_sample{$_};
	
	}

	#print OUTC join("\t",@common_samples),"\n";
	my $common_mutated_samples = scalar@common_samples;
	print OUTALL "$common_mutated_samples\t";

	#Print XX mutated samples
	foreach (keys %XX_mutated_samples){
	
		#print OUTE "$_\n";
	}
	my @XX_mutated_samples = keys %XX_mutated_samples;
	my $XX_mutated_samples = scalar@XX_mutated_samples;
	print OUTALL "$XX_mutated_samples\t";

	#Print interactor mutated samples
	foreach (keys %interactor_mutated_sample){
	
		#print OUTG "$_\n";
	}
	my @interactor_mutated_samples = keys %interactor_mutated_sample;
	my $interactor_mutated_sample = scalar@interactor_mutated_samples;
	print OUTALL "$interactor_mutated_sample\t"; 


	#Print samples with no G1 or G2 mutations 
	my @count;
	foreach (keys %unique_samples){
		
		$_=~s/\n|\t|\s//g;
		push (@count, $_) if exists $gene_sample_pairs{$gene_name}{$_};
		push (@count, $_) if exists $gene_sample_pairs{$interactor_name}{$_};
	
	}		
	my @mutated_unique_list = uniq @count;
	my $mutated_unique_list = scalar @mutated_unique_list;
	
	my $total_samples = keys %unique_samples;
	my $cases_altered_in_nethier_genes = $total_samples-$mutated_unique_list;
	print OUTALL "$cases_altered_in_nethier_genes\t";
	#print OUTALL scalar@mutated_unique_list,"\n";	
	
	
	#Calculate odds ratio

	my $cases_altered_in_XX_not_G2 = $XX_mutated_samples-$common_mutated_samples;
	my $cases_altered_in_G2_not_XX = $interactor_mutated_sample-$common_mutated_samples;
	
	next if (($common_mutated_samples == 0) || ($cases_altered_in_nethier_genes == 0) || ($cases_altered_in_XX_not_G2 == 0)|| ($cases_altered_in_G2_not_XX == 0));
			
	my $odds_ratio = (($common_mutated_samples* $cases_altered_in_nethier_genes)/($cases_altered_in_XX_not_G2*$cases_altered_in_G2_not_XX));
	my $odds_ratio_rounded =sprintf("%.3f",$odds_ratio);
	print OUTALL "$odds_ratio_rounded\t";
	my $log_odds_ratio_rounded = sprintf("%.3f",log10($odds_ratio));
	print OUTALL "$log_odds_ratio_rounded\t";

	#Calculate p-value on ODds ratio using Fischer test

	my $R = Statistics::R->new();
	$R->set( 'x', $common_mutated_samples);
	$R->set( 'y', $cases_altered_in_nethier_genes);
	$R->set( 'z', $cases_altered_in_XX_not_G2);
	$R->set( 'w', $cases_altered_in_G2_not_XX);
	$R->run('tab <-matrix(c(x,z,w,y), nrow=2)');
	$R->run( 'a <- fisher.test(tab)');
	my $squares = $R->get('a');
	#my @train = @$squares;
	my @fis = split("\n", $squares);
	$fis[3]=~s/p-value =//g;
	print OUTALL "$fis[3]\t$fis[9]";
	
	undef %XX_mutated_samples; undef %interactor_mutated_sample; undef %samples_not_altered_G1_G2; undef %unique_samples; undef %gene_sample_pairs;

}

#subroutine to calculate log10 of odds ratio
sub log10 {
	  my $n = shift;
	  return log($n)/log(10);
}

