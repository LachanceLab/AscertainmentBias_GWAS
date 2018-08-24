#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/shuffle/;
use Scalar::Util qw/looks_like_number/;

my $LOWER_TAIL_ONE = 7.5;
my $UPPER_TAIL_ZERO = 20;

my $start_run = time();

my $temp1; my @temp2=(); 
my @EAS_DAF=(); my @AMR_DAF=(); my @AFR_DAF=(); my @EUR_DAF=(); my @SAS_DAF=();
my @EAS_GRR=(); my @AMR_GRR=(); my @AFR_GRR=(); my @EUR_GRR=(); my @SAS_GRR=();
my @genotype_relative_risk_array = ();
#my $derived_freq; my $ancestral_freq;
my $i; my $s;
my $derived_power; my $ancestral_power;
my $derived_freq; my $ancestral_freq;
my $genotype_relative_risk; my $der_grr_power; my $anc_grr_power;

#****************************************************************************************************************#
# 								READING IN INPUT FILE TO RUN THE GAS POWER CALCULATOR							 #
# INPUT REQUIRES A FORMAT OF:																					 #
# Data: PATH TO DATA FILE																						 #
# Number of simulations: #																						 #
# Number of randomly selected SNPs: #																			 #
# Study population: STUDY POPULATIONS																			 #
# Genotype relative risk: #																						 #
# Prevalence: # 																						 		 #
# Number of cases: #																						 	 #
# Number of controls: #																						 	 #
# P-value: #																						 			 #
# Model: ADDITIVE/MULTIPLICATIVE/RECESSIVE/DOMINANT																 #
#****************************************************************************************************************#

my $input_temp1 = `cat gas_power_inputs.txt`;
#print $input_temp, "\n";


my @input_temp2 = split /\n/, $input_temp1;

my $data; my $number_of_simulations; my $num_randomly_selected_SNPs; my $study_pop; 
my $genotype_relative_risk_temp; my $prev; my $ncases; my $ncontrols; my $alpha; my $dominance;

for ($i=0; $i<scalar(@input_temp2); $i++) {
	my @input_temp3 = split /: /, $input_temp2[$i];
	#print $input_temp3[1], "\n";
	
	if ($i == 0) {
		$data = $input_temp3[1];
	}
	elsif ($i == 1) {
		$number_of_simulations = $input_temp3[1];
	}
	elsif ($i == 2) {
		$num_randomly_selected_SNPs = $input_temp3[1];
	}
	elsif ($i == 3) {
		$study_pop = $input_temp3[1];
	}
	elsif ($i == 4) {
		$genotype_relative_risk_temp = $input_temp3[1];
	}
	elsif ($i == 5) {
		$prev = $input_temp3[1];
	}
	elsif ($i == 6) {
		$ncases = $input_temp3[1];
	}
	elsif ($i == 7) {
		$ncontrols = $input_temp3[1];
	}
	elsif ($i == 8) {
		$alpha = $input_temp3[1];
	}
	elsif ($i == 9) {
		$dominance = $input_temp3[1];
	}
}


if ($study_pop eq "eas" || $study_pop eq "EAS" || $study_pop eq "Eas") {
	$temp1 = `tail -n +2 $data | awk '\$7!=0' | awk '\$7!=1'`;
}
elsif ($study_pop eq "amr" || $study_pop eq "AMR" || $study_pop eq "Amr") {
	$temp1 = `tail -n +2 $data | awk '\$8!=0' | awk '\$8!=1'`;
}
elsif ($study_pop eq "afr" || $study_pop eq "AFR" || $study_pop eq "Afr") {
	$temp1 = `tail -n +2 $data | awk '\$9!=0' | awk '\$9!=1'`;
}
elsif ($study_pop eq "eur" || $study_pop eq "EUR" || $study_pop eq "Eur") {
	$temp1 = `tail -n +2 $data | awk '\$10!=0' | awk '\$10!=1'`;
}
elsif ($study_pop eq "sas" || $study_pop eq "SAS" || $study_pop eq "Sas") {
	$temp1 = `tail -n +2 $data | awk '\$11!=0' | awk '\$11!=1'`;
}
#print $temp1, "\n";
my $midpoint_run = time();
my $mid_run_time = $midpoint_run - $start_run;
#print "mid run time = ", $mid_run_time, "\n";

@temp2 = split /\n/, $temp1;
#print $temp2[0], "\n";

my $sim;
open my $fh_output, '>', "Simulation_output.txt" or die;
print $fh_output "Rep", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";


for ($sim=0; $sim<$number_of_simulations; $sim++) {
	
	`Rscript gamma_distributed_OR_generator.R`;					#CALLING RSCRIPT FOR GENERATING GENETIC RELATIVE RISK (GRR)
	my $GRR_file = "AFR_EUR_GRR.txt";							#AFR_EUR_GRR.TXT IS THE NAME OF THE FILE THAT RSCRIPT CREATED
	
	my $GRR_file_nrows;
	open (GRR_fh_in, "<". $GRR_file);
	$GRR_file_nrows++ while (<GRR_fh_in>);
	close (GRR_fh_in);
	
	
	my $GRR; my @GRR_temp = (); my @GRR_list = (); my @other_GRR = ();
	open (GRR_fh_in, "<". $GRR_file);
	for ($s=0; $s<$GRR_file_nrows; $s++) {
		$GRR = <GRR_fh_in>;
		@GRR_temp = split("\t", $GRR);
		$GRR_list[$s] = $GRR_temp[0];
		$other_GRR[$s] = $GRR_temp[1];
		
		#print $GRR_list[$s], "\n";
		#print $other_GRR[$s], "\n";
	}
	close (GRR_fh_in);
	
	#print @GRR_list, "\n";

	my @power_derived=(); my @power_ancestral=();
	my @der_power_grr=(); my @anc_power_grr=();
	
	open my $fh_anc_temp, '>', "anc_temp.txt" or die;
	print $fh_anc_temp "power", "\t", "AFR_AAF", "\t", "AMR_AAF", "\t", "EAS_AAF", "\t", "EUR_AAF", "\t", "SAS_AAF", "\t", "AFR_GRR", "\t", "EUR_GRR", "\n";
	
	open my $fh_der_temp, '>', "der_temp.txt" or die;
	print $fh_der_temp "power", "\t", "AFR_DAF", "\t", "AMR_DAF", "\t", "EAS_DAF", "\t", "EUR_DAF", "\t", "SAS_DAF", "\t", "AFR_GRR", "\t", "EUR_GRR", "\n";
	
	
	for ($i=0; $i<scalar(@temp2); $i++) {
		
		my @temp3=();
		@temp3 = split /\t/, $temp2[$i];

		$EAS_DAF[$i] = $temp3[6];
		$AMR_DAF[$i] = $temp3[7];
		$AFR_DAF[$i] = $temp3[8];
		$EUR_DAF[$i] = $temp3[9];
		$SAS_DAF[$i] = $temp3[10];
		$genotype_relative_risk = $GRR_list[$i];
		
		#print $genotype_relative_risk, "\n";
		
		
		$EAS_GRR[$i] = $other_GRR[$i];
		$AMR_GRR[$i] = $other_GRR[$i];
		$AFR_GRR[$i] = $GRR_list[$i];
		$EUR_GRR[$i] = $other_GRR[$i];
		$SAS_GRR[$i] = $other_GRR[$i];
		
		

		if ($study_pop eq "eas" || $study_pop eq "EAS" || $study_pop eq "Eas") {
			$derived_freq = $EAS_DAF[$i];
			$ancestral_freq = 1 - $EAS_DAF[$i];
		}
		elsif ($study_pop eq "amr" || $study_pop eq "AMR" || $study_pop eq "Amr") {
			$derived_freq = $AMR_DAF[$i];
			$ancestral_freq = 1 - $AMR_DAF[$i];
		}
		elsif ($study_pop eq "afr" || $study_pop eq "AFR" || $study_pop eq "Afr") {
			$derived_freq = $AFR_DAF[$i];
			$ancestral_freq = 1 - $AFR_DAF[$i];
		}
		elsif ($study_pop eq "eur" || $study_pop eq "EUR" || $study_pop eq "Eur") {
			$derived_freq = $EUR_DAF[$i];
			$ancestral_freq = 1 - $EUR_DAF[$i];
		}
		elsif ($study_pop eq "sas" || $study_pop eq "SAS" || $study_pop eq "Sas") {
			$derived_freq = $SAS_DAF[$i];
			$ancestral_freq = 1 - $SAS_DAF[$i];
		}
		

		derived_calculate ($ncases, $ncontrols, $derived_freq, $genotype_relative_risk, $prev, $alpha);				#CALLING GAS POWER CALCULATOR
		$power_derived[$i] = $derived_power - 1;
		
		ancestral_calculate ($ncases, $ncontrols, $ancestral_freq, $genotype_relative_risk, $prev, $alpha);			#CALLING GAS POWER CALCULATOR
		$power_ancestral[$i] = $ancestral_power - 1;
		
		der_grr_power_calculate ($ncases, $ncontrols, $derived_freq, $genotype_relative_risk, $prev, $alpha);		#CALLING GAS POWER CALCULATOR
		$der_power_grr[$i] = $der_grr_power - 1;
		
		anc_grr_power_calculate ($ncases, $ncontrols, $ancestral_freq, $genotype_relative_risk, $prev, $alpha);		#CALLING GAS POWER CALCULATOR
		$anc_power_grr[$i] = $anc_grr_power - 1;
	
		#print "power derived=", $power_derived[$i], "\n";
		#print "power ancestral=", $power_ancestral[$i], "\n";
		#print "power grr=", $power_grr[$i], "\n";
		
		print $fh_anc_temp $anc_power_grr[$i], "\t", (1-$AFR_DAF[$i]), "\t", (1-$AMR_DAF[$i]), "\t", (1-$EAS_DAF[$i]), "\t", (1-$EUR_DAF[$i]), "\t", (1-$SAS_DAF[$i]), "\t", $GRR_list[$i], "\t", $other_GRR[$i];
		
		
		print $fh_der_temp $der_power_grr[$i], "\t", $AFR_DAF[$i], "\t", $AMR_DAF[$i], "\t", $EAS_DAF[$i], "\t", $EUR_DAF[$i], "\t", $SAS_DAF[$i], "\t", $GRR_list[$i], "\t", $other_GRR[$i];


	}
	close ($fh_anc_temp);
	close ($fh_der_temp);
	
	my $data2; my $data3;
	$data2 = "anc_temp.txt";
	$data3 = "der_temp.txt";
	
	my $temp4; my $temp5;
	
	$temp4 = `tail -n +2 $data2 | awk '\$1!=0'`;			#SELECTING SNPS THAT WERE DETECTED BY THE POWER CALCULATOR
	$temp5 = `tail -n +2 $data3 | awk '\$1!=0'`;
	
	my @temp6 = (); my @temp7 = ();
	@temp6 = split /\n/, $temp4;
	@temp7 = split /\n/, $temp5;
	
	
	my @shuffled_temp6=(); my @shuffled_temp7=();
	@shuffled_temp6 = shuffle(@temp6);
	@shuffled_temp7 = shuffle(@temp7);
	
	my $power_derived_sum=0; my $power_ancestral_sum=0;
	my $der_power_grr_sum=0; my $anc_power_grr_sum=0;
	
	
	my $derived_weighted_count_sum=0;
	my $derived_AFR_weighted_sum=0; my @derived_AFR_weighted=();
	my $derived_AMR_weighted_sum=0; my @derived_AMR_weighted=();
	my $derived_EAS_weighted_sum=0; my @derived_EAS_weighted=();
	my $derived_EUR_weighted_sum=0; my @derived_EUR_weighted=();
	my $derived_SAS_weighted_sum=0; my @derived_SAS_weighted=();
	
	my $ancestral_weighted_count_sum=0;
	my $ancestral_AFR_weighted_sum=0; my @ancestral_AFR_weighted=();
	my $ancestral_AMR_weighted_sum=0; my @ancestral_AMR_weighted=();
	my $ancestral_EAS_weighted_sum=0; my @ancestral_EAS_weighted=();
	my $ancestral_EUR_weighted_sum=0; my @ancestral_EUR_weighted=();
	my $ancestral_SAS_weighted_sum=0; my @ancestral_SAS_weighted=();

	my @ancestral_post_gwas_AFR_GRR = (); my @ancestral_post_gwas_EUR_GRR = ();
	my @derived_post_gwas_AFR_GRR =(); my @derived_post_gwas_EUR_GRR =();
	my @ancestral_power = (); my @derived_power = ();
	my @AFR_AAF = (); my @AMR_AAF = (); my @EAS_AAF = (); my @EUR_AAF = (); my @SAS_AAF = ();
	my @AFR_DAF1 = (); my @AMR_DAF1 = (); my @EAS_DAF1 = (); my @EUR_DAF1 = (); my @SAS_DAF1 = ();
	
	
	open my $anc_post_gwas_grr_out, '>', "GRR_post_GWAS_ancestral.txt" or die;
	print $anc_post_gwas_grr_out "AFR_GRR_post_GWAS", "\t", "EUR_GRR_post_GWAS", "\n";
	
	open my $der_post_gwas_grr_out, '>', "GRR_post_GWAS_derived.txt" or die;
	print $der_post_gwas_grr_out "AFR_GRR_post_GWAS", "\t", "EUR_GRR_post_GAS", "\n";
	
	for ($i=0; $i<$num_randomly_selected_SNPs; $i++) {
		my @temp8=(); my @temp9=();
		@temp8 = split /\t/, $shuffled_temp6[$i];
		@temp9 = split /\t/, $shuffled_temp7[$i];

		$ancestral_power[$i] = $temp8[0];
		$EAS_AAF[$i] = $temp8[1];
		$AMR_AAF[$i] = $temp8[2];
		$AFR_AAF[$i] = $temp8[3];
		$EUR_AAF[$i] = $temp8[4];
		$SAS_AAF[$i] = $temp8[5];
		$ancestral_post_gwas_AFR_GRR[$i] = $temp8[6];
		$ancestral_post_gwas_EUR_GRR[$i] = $temp8[7];
		
		$derived_power[$i] = $temp9[0];
		$AFR_DAF1[$i] = $temp9[1];
		$AMR_DAF1[$i] = $temp9[2];
		$EAS_DAF1[$i] = $temp9[3];
		$EUR_DAF1[$i] = $temp9[4];
		$SAS_DAF1[$i] = $temp9[5];
		$derived_post_gwas_AFR_GRR[$i] = $temp9[6];
		$derived_post_gwas_EUR_GRR[$i] = $temp9[7];
		
		$power_derived_sum = $power_derived_sum + $derived_power[$i];
		$power_ancestral_sum = $power_ancestral_sum + $ancestral_power[$i];
		
		
	
	
		my $derived_weighted_count = $derived_power[$i] / $power_derived_sum;
		$derived_weighted_count_sum = $derived_weighted_count_sum + $derived_weighted_count;
		
		$derived_AFR_weighted[$i] = $derived_weighted_count * $AFR_DAF1[$i];
		$derived_AFR_weighted_sum = $derived_AFR_weighted_sum + $derived_AFR_weighted[$i];
		
		$derived_AMR_weighted[$i] = $derived_weighted_count * $AMR_DAF1[$i];
		$derived_AMR_weighted_sum = $derived_AMR_weighted_sum + $derived_AMR_weighted[$i];
		
		$derived_EAS_weighted[$i] = $derived_weighted_count * $EAS_DAF1[$i];
		$derived_EAS_weighted_sum = $derived_EAS_weighted_sum + $derived_EAS_weighted[$i];

		$derived_EUR_weighted[$i] = $derived_weighted_count * $EUR_DAF1[$i];
		$derived_EUR_weighted_sum = $derived_EUR_weighted_sum + $derived_EUR_weighted[$i];
		
		$derived_SAS_weighted[$i] = $derived_weighted_count * $SAS_DAF1[$i];
		$derived_SAS_weighted_sum = $derived_SAS_weighted_sum + $derived_SAS_weighted[$i];

		#print "Power = ", $power_derived[$i], "\n";
		#print "EAS DAF = ", $EAS_DAF[$i], "\n";
		
		my $ancestral_weighted_count = $ancestral_power[$i] / $power_ancestral_sum;
		$ancestral_weighted_count_sum = $ancestral_weighted_count_sum + $ancestral_weighted_count;
		
		$ancestral_AFR_weighted[$i] = $ancestral_weighted_count * $AFR_AAF[$i];
		$ancestral_AFR_weighted_sum = $ancestral_AFR_weighted_sum + $ancestral_AFR_weighted[$i];
		
		$ancestral_AMR_weighted[$i] = $ancestral_weighted_count * $AMR_AAF[$i];
		$ancestral_AMR_weighted_sum = $ancestral_AMR_weighted_sum + $ancestral_AMR_weighted[$i];
		
		$ancestral_EAS_weighted[$i] = $ancestral_weighted_count * $EAS_AAF[$i];
		$ancestral_EAS_weighted_sum = $ancestral_EAS_weighted_sum + $ancestral_EAS_weighted[$i];
		
		$ancestral_EUR_weighted[$i] = $ancestral_weighted_count * $EUR_AAF[$i];
		$ancestral_EUR_weighted_sum = $ancestral_EUR_weighted_sum + $ancestral_EUR_weighted[$i];
		
		$ancestral_SAS_weighted[$i] = $ancestral_weighted_count * $SAS_AAF[$i];
		$ancestral_SAS_weighted_sum = $ancestral_SAS_weighted_sum + $ancestral_SAS_weighted[$i];
	
	
		print $anc_post_gwas_grr_out $ancestral_post_gwas_AFR_GRR[$i], "\t", $ancestral_post_gwas_EUR_GRR[$i], "\n";
		print $der_post_gwas_grr_out $derived_post_gwas_AFR_GRR[$i], "\t", $derived_post_gwas_AFR_GRR[$i], "\n";
	}
	close ($anc_post_gwas_grr_out);
	close ($der_post_gwas_grr_out);
	
	my $AFR_DER = $derived_AFR_weighted_sum / $derived_weighted_count_sum;
	my $AMR_DER = $derived_AMR_weighted_sum / $derived_weighted_count_sum;
	my $EAS_DER = $derived_EAS_weighted_sum / $derived_weighted_count_sum;
	my $EUR_DER = $derived_EUR_weighted_sum / $derived_weighted_count_sum;
	my $SAS_DER = $derived_SAS_weighted_sum / $derived_weighted_count_sum;
	
	my $AFR_ANC = $ancestral_AFR_weighted_sum / $ancestral_weighted_count_sum;
	my $AMR_ANC = $ancestral_AMR_weighted_sum / $ancestral_weighted_count_sum;
	my $EAS_ANC = $ancestral_EAS_weighted_sum / $ancestral_weighted_count_sum;
	my $EUR_ANC = $ancestral_EUR_weighted_sum / $ancestral_weighted_count_sum;
	my $SAS_ANC = $ancestral_SAS_weighted_sum / $ancestral_weighted_count_sum;

	print $fh_output $sim+1, "\t", $AFR_ANC, "\t", $AMR_ANC, "\t", $EAS_ANC, "\t", $EUR_ANC, "\t", $SAS_ANC, "\t", $AFR_DER, "\t", $AMR_DER, "\t", $EAS_DER, "\t", $EUR_DER, "\t", $SAS_DER, "\n";

	print $sim+1, "\n";

}
close ($fh_output);

	
	
	



#***************************************************#
# GAS Power (from University of Michigan) Calculator#
# From GAS-power-calculator GitHub, INC.			#
#***************************************************#


sub ninv {

	my $split1; my $split2; my $const1; my $const2; 
	my @aa=(); my @bb=(); my @cc=(); my @dd=(); my @ee=(); my @ff=();
	my $prob; my $q; my $r; my $x;


	$split1 = 0.425;
	$split2 = 5.0;
	$const1 = 0.180625;
	$const2 = 1.6;
	
	@aa=(3.3871328727963666080, 133.14166789178437745, 1971.5909503065514427, 13731.693765509461125,
      45921.953931549871457, 67265.770927008700853, 33430.575583588128105, 2509.0809287301226727);
      
	@bb=(42.313330701600911252, 687.18700749205790830, 5394.1960214247511077, 21213.794301586595867, 
      39307.895800092710610, 28729.085735721942674, 5226.4952788528545610);
    
	@cc=(1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504,
      1.27045825245236838258, 0.241780725177450611770, 0.0227238449892691845833, 0.000774545014278341407640);
    
	@dd=(2.05319162663775882187, 1.67638483018380384940, 0.689767334985100004550, 0.148103976427480074590, 
      0.0151986665636164571966, 0.000547593808499534494600, 0.00000000105075007164441684324);
      
	@ee=(6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 0.296560571828504891230,
      0.0265321895265761230930, 0.00124266094738807843860, 0.0000271155556874348757815, 0.000000201033439929228813265);
    
	@ff=(0.599832206555887937690, 0.136929880922735805310, 0.0148753612908506148525, 0.000786869131145613259100, 
      0.0000184631831751005468180, 0.000000142151175831644588870, 0.00000000000000204426310338993978564);
    
    
    $prob = $alpha * 0.5;
    $q = $prob - 0.5;
    
    if (abs($q) < $split1) {
    	$r = $const1 - $q * $q;
    	return $q * ((((((( $aa[7] * $r + $aa[6] ) * $r + $aa[5] ) * $r + $aa[4] ) * $r + $aa[3] ) * $r + $aa[2] ) * $r + $aa[1] ) * $r + $aa[0] ) /
    				((((((( $bb[6] * $r + $bb[5] ) * $r + $bb[4] ) * $r + $bb[3] ) * $r + $bb[2] ) * $r + $bb[1] ) * $r + $bb[0] ) * $r * 1.0 );
    }
    else {
    	if ($q < 0) {
    		$r = $prob;
    	}
    	else {
    		$r = 1.0 - $prob;
    	}
    	if ($r < 0.0000000001) {
    		return ($q < 0 ? -20.0 : 20.0);
    	}
    	if ($r > 0.0) {
    		$r = sqrt( -1 * log($r) );
    		if ($r <= $split2) {
    			$r = $r - $const2;
    			$x = ((((((( $cc[7] * $r + $cc[6] ) * $r + $cc[5] ) * $r + $cc[4] ) * $r + $cc[3] ) * $r + $cc[2] ) * $r + $cc[1] ) * $r + $cc[0] ) /
    				 ((((((( $dd[6] * $r + $dd[5] ) * $r + $dd[4] ) * $r + $dd[3] ) * $r + $dd[2] ) * $r + $dd[1] ) * $r + $dd[0] ) * $r + 1 );
    		}
    		else {
    			$r = $r - $split2;
    			$x = ((((((( $ee[7] * $r + $ee[6] ) * $r + $ee[5] ) * $r + $ee[4] ) * $r + $ee[3] ) * $r + $ee[2] ) * $r + $ee[1] ) * $r + $ee[0] ) /
    				 ((((((( $ff[6] * $r + $ff[5] ) * $r + $ff[4] ) * $r + $ff[3] ) * $r + $ff[2] ) * $r + $ff[1] ) * $r + $ff[0] ) * $r + 1 );
    		}
    	}
    	else {
    		$x = 9;
    	}
    	if ($q < 0) {
    		$x = -1 * $x;
    	}
    	return $x;
    }
}


my $upper; my $p;

sub ndist {

	my $z; my $y; 

	$z = shift;
	$upper = shift;
	
	if ($z < 0) {
		$upper = 'false';
		$z = -1 * $z;
	}
	
	if ( (($z > $LOWER_TAIL_ONE) && ($upper eq 'false')) || ($z > $UPPER_TAIL_ZERO) ) {
		return ($upper eq 'true' ? 0.0 : 1.0);
	}
	
	$y = 0.5 * $z * $z;
	
	if ($z < 1.28) {
		$p = 0.5 - $z * (0.398942280444 - 0.399903438504 * $y /
			($y + 5.75885480458 - 29.8213557808 /
			($y + 2.62433121679 + 48.6959930692 /
			($y + 5.92885724438))));
	}
	else {
		$p =  0.398942270385 * exp($y * -1) /
			($z - 0.000000028052 + 1.00000615302 /
			($z + 0.0003980648 + 1.98615381364 /
			($z - 0.151679116635 + 5.29330324926 /
			($z + 4.8385912808 - 15.1508972451 /
			($z + 0.742380924027 + 30.789933034 /
			($z + 3.99019417011))))));
	}
	return ($upper eq 'true' ? $p : 1 - $p);
}


sub derived_calculate {

	my @pp=(); my $aa_freq; my $ab_freq; my $bb_freq; my @f=();
	my @additive=(); my @multiplicative=(); my @dominant=(); my @recessive=();
	my $scale; my $aa_prob; my $ab_prob; my $bb_prob;
	my $C; my $pcases; my $pcontrols; my $vcases; my $vcontrols; my $ncp;


	
	@pp = ($derived_freq * $derived_freq, 2 * $derived_freq * (1 - $derived_freq), (1. - $derived_freq) * (1. - $derived_freq));

	$aa_freq = $pp[0];
	$ab_freq = $pp[1];
	$bb_freq = $pp[2];

	if ($dominance eq "additive") {
		@additive = (2 * $genotype_relative_risk - 1.0, $genotype_relative_risk, 1.0);
		@f = @additive;
	}
	elsif ($dominance eq "multiplicative") {
		@multiplicative = ($genotype_relative_risk * $genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @multiplicative;
	}
	elsif ($dominance eq "dominant") {
		@dominant = ($genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @dominant;
	}
	elsif ($dominance eq "recessive") {
		@recessive = ($genotype_relative_risk, 1.0, 1.0);
		@f = @recessive;
	}
	
	$scale = $prev / ($f[0] * $aa_freq + $f[1] * $ab_freq + $f[2] * $bb_freq);
	$f[0] = $f[0] * $scale;
	$f[1] = $f[1] * $scale;
	$f[2] = $f[2] * $scale;

	$aa_prob = $f[0];
	$ab_prob = $f[1];
	$bb_prob = $f[2];

	$C = -1 * ninv();
	$pcases = ($f[0] * $pp[0] + $f[1] * $pp[1] * 0.5) / $prev;
	$pcontrols = ( (1. - $f[0]) * $pp[0] + (1. - $f[1]) * $pp[1] * 0.5) / (1. - $prev);
	$vcases = $pcases * (1.0 - $pcases);
	$vcontrols = $pcontrols * (1.0 - $pcontrols);
	$ncp = ($pcases - $pcontrols) / sqrt( ($vcases / $ncases + $vcontrols / $ncontrols) * 0.5 );
	$derived_power = ndist( -1 * $C - $ncp, 'false') + ndist($C - $ncp, 'true');
	return $derived_power;
}


sub ancestral_calculate {

	my @pp=(); my $aa_freq; my $ab_freq; my $bb_freq; my @f=();
	my @additive=(); my @multiplicative=(); my @dominant=(); my @recessive=();
	my $scale; my $aa_prob; my $ab_prob; my $bb_prob;
	my $C; my $pcases; my $pcontrols; my $vcases; my $vcontrols; my $ncp;


	
	@pp = ($ancestral_freq * $ancestral_freq, 2 * $ancestral_freq * (1 - $ancestral_freq), (1. - $ancestral_freq) * (1. - $ancestral_freq));

	$aa_freq = $pp[0];
	$ab_freq = $pp[1];
	$bb_freq = $pp[2];

	if ($dominance eq "additive") {
		@additive = (2 * $genotype_relative_risk - 1.0, $genotype_relative_risk, 1.0);
		@f = @additive;
	}
	elsif ($dominance eq "multiplicative") {
		@multiplicative = ($genotype_relative_risk * $genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @multiplicative;
	}
	elsif ($dominance eq "dominant") {
		@dominant = ($genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @dominant;
	}
	elsif ($dominance eq "recessive") {
		@recessive = ($genotype_relative_risk, 1.0, 1.0);
		@f = @recessive;
	}
	
	$scale = $prev / ($f[0] * $aa_freq + $f[1] * $ab_freq + $f[2] * $bb_freq);
	$f[0] = $f[0] * $scale;
	$f[1] = $f[1] * $scale;
	$f[2] = $f[2] * $scale;

	$aa_prob = $f[0];
	$ab_prob = $f[1];
	$bb_prob = $f[2];

	$C = -1 * ninv();
	$pcases = ($f[0] * $pp[0] + $f[1] * $pp[1] * 0.5) / $prev;
	$pcontrols = ( (1. - $f[0]) * $pp[0] + (1. - $f[1]) * $pp[1] * 0.5) / (1. - $prev);
	$vcases = $pcases * (1.0 - $pcases);
	$vcontrols = $pcontrols * (1.0 - $pcontrols);
	$ncp = ($pcases - $pcontrols) / sqrt( ($vcases / $ncases + $vcontrols / $ncontrols) * 0.5 );
	$ancestral_power = ndist( -1 * $C - $ncp, 'false') + ndist($C - $ncp, 'true');
	return $ancestral_power;
}



sub der_grr_power_calculate {
	my @pp=(); my $aa_freq; my $ab_freq; my $bb_freq; my @f=();
	my @additive=(); my @multiplicative=(); my @dominant=(); my @recessive=();
	my $scale; my $aa_prob; my $ab_prob; my $bb_prob;
	my $C; my $pcases; my $pcontrols; my $vcases; my $vcontrols; my $ncp;


	
	@pp = ($derived_freq * $derived_freq, 2 * $derived_freq * (1 - $derived_freq), (1. - $derived_freq) * (1. - $derived_freq));

	$aa_freq = $pp[0];
	$ab_freq = $pp[1];
	$bb_freq = $pp[2];

	if ($dominance eq "additive") {
		@additive = (2 * $genotype_relative_risk - 1.0, $genotype_relative_risk, 1.0);
		@f = @additive;
	}
	elsif ($dominance eq "multiplicative") {
		@multiplicative = ($genotype_relative_risk * $genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @multiplicative;
	}
	elsif ($dominance eq "dominant") {
		@dominant = ($genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @dominant;
	}
	elsif ($dominance eq "recessive") {
		@recessive = ($genotype_relative_risk, 1.0, 1.0);
		@f = @recessive;
	}
	
	$scale = $prev / ($f[0] * $aa_freq + $f[1] * $ab_freq + $f[2] * $bb_freq);
	$f[0] = $f[0] * $scale;
	$f[1] = $f[1] * $scale;
	$f[2] = $f[2] * $scale;

	$aa_prob = $f[0];
	$ab_prob = $f[1];
	$bb_prob = $f[2];

	$C = -1 * ninv();
	$pcases = ($f[0] * $pp[0] + $f[1] * $pp[1] * 0.5) / $prev;
	$pcontrols = ( (1. - $f[0]) * $pp[0] + (1. - $f[1]) * $pp[1] * 0.5) / (1. - $prev);
	$vcases = $pcases * (1.0 - $pcases);
	$vcontrols = $pcontrols * (1.0 - $pcontrols);
	$ncp = ($pcases - $pcontrols) / sqrt( ($vcases / $ncases + $vcontrols / $ncontrols) * 0.5 );
	$der_grr_power = ndist( -1 * $C - $ncp, 'false') + ndist($C - $ncp, 'true');
	return $der_grr_power;
}


sub anc_grr_power_calculate {
	my @pp=(); my $aa_freq; my $ab_freq; my $bb_freq; my @f=();
	my @additive=(); my @multiplicative=(); my @dominant=(); my @recessive=();
	my $scale; my $aa_prob; my $ab_prob; my $bb_prob;
	my $C; my $pcases; my $pcontrols; my $vcases; my $vcontrols; my $ncp;


	
	@pp = ($ancestral_freq * $ancestral_freq, 2 * $ancestral_freq * (1 - $ancestral_freq), (1. - $ancestral_freq) * (1. - $ancestral_freq));

	$aa_freq = $pp[0];
	$ab_freq = $pp[1];
	$bb_freq = $pp[2];

	if ($dominance eq "additive") {
		@additive = (2 * $genotype_relative_risk - 1.0, $genotype_relative_risk, 1.0);
		@f = @additive;
	}
	elsif ($dominance eq "multiplicative") {
		@multiplicative = ($genotype_relative_risk * $genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @multiplicative;
	}
	elsif ($dominance eq "dominant") {
		@dominant = ($genotype_relative_risk, $genotype_relative_risk, 1.0);
		@f = @dominant;
	}
	elsif ($dominance eq "recessive") {
		@recessive = ($genotype_relative_risk, 1.0, 1.0);
		@f = @recessive;
	}
	
	$scale = $prev / ($f[0] * $aa_freq + $f[1] * $ab_freq + $f[2] * $bb_freq);
	$f[0] = $f[0] * $scale;
	$f[1] = $f[1] * $scale;
	$f[2] = $f[2] * $scale;

	$aa_prob = $f[0];
	$ab_prob = $f[1];
	$bb_prob = $f[2];

	$C = -1 * ninv();
	$pcases = ($f[0] * $pp[0] + $f[1] * $pp[1] * 0.5) / $prev;
	$pcontrols = ( (1. - $f[0]) * $pp[0] + (1. - $f[1]) * $pp[1] * 0.5) / (1. - $prev);
	$vcases = $pcases * (1.0 - $pcases);
	$vcontrols = $pcontrols * (1.0 - $pcontrols);
	$ncp = ($pcases - $pcontrols) / sqrt( ($vcases / $ncases + $vcontrols / $ncontrols) * 0.5 );
	$anc_grr_power = ndist( -1 * $C - $ncp, 'false') + ndist($C - $ncp, 'true');
	return $anc_grr_power;
}	
	
	
	
	
	