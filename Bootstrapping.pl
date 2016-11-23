# perform t.test of wild type and mutant samples. Input data (gene expression matrix) should be sorted according to mutation information (wildtype to mutant as 0/1)
# Keeps 10% of wild and mutant samples from dataset and perform bootstrapping to get difference in mean expression values
# Swati Kaushik 2016-02-16

use strict;
use warnings;
use Statistics::R;

my $filename = $ARGV[0];
open GE, "$filename" or die $!;

my @exprfile =<GE>;
shift @exprfile;

#write output
open OUT, ">>muliple-test-bootstrapping" or die $!;

foreach (@exprfile){

	my @rows = split("\t", $_);
	my $geneid = shift @rows;
	my $R = Statistics::R->new();
	$R->set( 'x', \@rows);
	$R->run( 'a <- t.test(x[2:527],x[528:577])[c("p.value","estimate")]');
	my $squaresn = $R->get('a');
	my @trainn = @$squaresn;
	print OUT "$geneid\t$trainn[2]\t$trainn[10]\t$trainn[11]\t";
	

	my $n=0; my @train_pvalue; my @test_pvalue; my @train_w_mean; my @train_m_mean;  my @test_w_mean; my @test_m_mean;

	while ($n <100)
	{
		#random number generation for wild type samples 
		my $maxwild = 525; 	my $minwild = 0; my @randomwild; my @randomwild_GE_training; my %wildhash;my @randomwild_GE_test; 
		
		for (my $i=0; $i<53; $i++){
			
			#my $random_number_wild = int(rand($maxwild)) + $minwild;
			my $random_number_wild = int($minwild + rand($maxwild - $minwild));
			$wildhash{$random_number_wild}=$rows[$random_number_wild];
			#print OUT "$rows[$random_number_wild]\t";
			push (@randomwild_GE_test, $rows[$random_number_wild]);
		
		}
		
		for (my $i=0; $i<525; $i++){
			
			if (not exists ($wildhash{$i}))
				{	push (@randomwild_GE_training, $rows[$i]);
					#print OUT "$i\t$rows[$i]\n";
				}
		
		}
		
		#random number generation for mutant samples 
		my $maxmutant = 575; my $minmutant = 526;	my @randommutant; my @randommutant_GE_training; my %mutanthash; my@randommutant_GE_test;
		
		for (my $i=0; $i<5; $i++){
			
			my $random_number_mutant = int($minmutant + rand($maxmutant - $minmutant));
			#print "$random_number_mutant\n";
			if (!(exists($mutanthash{$random_number_mutant})))
				{
					$mutanthash{$random_number_mutant}=$rows[$random_number_mutant];
					push (@randommutant_GE_test, $rows[$random_number_mutant]);
					#print OUT "**$rows[$random_number_mutant]\t";
				}
		}		
		
		for (my $i=526; $i<575; $i++){
			
			if (not exists ($mutanthash{$i}))
				{	push (@randommutant_GE_training, $rows[$i]);
					#print OUT "$i\t$rows[$i]\n";
				}
		
		}

	# t.test of training set using R

	my $R = Statistics::R->new();
	$R->set( 'x', \@randomwild_GE_training);
	$R->set( 'y', \@randommutant_GE_training);
	$R->run( 'a <- t.test(x,y)[c("p.value","estimate")]');
	my $squares = $R->get('a');
	my @train = @$squares;
	#print "$train[2]\t$train[10]\t$train[11]\t";
	push (@train_pvalue, $train[2]); push (@train_w_mean, $train[10]); push (@train_m_mean, $train[11]);
	
	# t.test of test set using R

	$R->set( 'x1', \@randomwild_GE_test);
	$R->set( 'y1', \@randommutant_GE_test);
	$R->run( 'a1 <- t.test(x1,y1)[c("p.value","estimate")]');
	my $squares1 = $R->get('a1');
	my @train1 = @$squares1;
	#print "$train1[2]\t$train1[10]\t$train1[11]\n";
	push (@test_pvalue, $train1[2]); push (@test_w_mean, $train1[10]); push (@test_m_mean, $train1[11]);

	$n++;	
	}
	
print OUT average(@train_pvalue), "\t";
print OUT average(@train_w_mean), "\t";
print OUT average(@train_m_mean), "\t";
print OUT average(@test_pvalue), "\t";
print OUT average(@test_w_mean), "\t";
print OUT average(@test_m_mean), "\n";

}

#sub routine to calculate average of pvalue and means

sub average {

	my @array =@_;
	my $sum;
	foreach (@array) {$sum+=$_;}
	return $sum/@array;
}
