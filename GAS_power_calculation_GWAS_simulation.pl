#***************************#
# Perl code by Michelle Kim #
#    mkim634@gatech.edu     #
#***************************#

#!/usr/bin/perl -w
use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use File::chdir;
use Getopt::Long;
use List::Util qw/shuffle/;

my $LOWER_TAIL_ONE = 7.5;
my $UPPER_TAIL_ZERO = 20;


#**********#
#GetOptions#
#**********#
my $working_directory; my $replicates; my $technology; my $cases; my $controls;
my $cases_file; my $controls_file;
my $study_pop; my $significance; my $prevalence;
my $significance_file; my $prevalence_file;
my $grr; my $dominance;
my $grr_file;
my $help;


GetOptions ('w=s' => \$working_directory,
			'n=i' => \$replicates,
			't=s' => \$technology,
			'cases=i' => \$cases,
			'casesfh=s' => \$cases_file,
			'controls=i' => \$controls,
			'controlsfh=s' => \$controls_file,
			'pop=s' => \$study_pop,
			's=f' => \$significance,
			'sfh=s' => \$significance_file,
			'p=f' => \$prevalence,
			'pfh=s' => \$prevalence_file,
			'r=f' => \$grr,
			'rfh=s' => \$grr_file,
			'd=s' => \$dominance,
			'help' => \$help
) or die("Error in command line arguments\n");

if ($help) {
	print "Usage: perl gwas.pl [options] \n";
	print "-w: working directory \n";
	print "-n: number of replicates \n";
	print "-t: technology used \n";
	print "-cases: number of cases \n";
	print "-casesfh: number of cases in file \n";
	print "-controls: number of controls \n";
	print "-controlsfh: number of controls file \n";
	print "-pop: study population \n";
	print "-s: significance (alpha) \n";
	print "-p: prevalence \n";
	print "-r: genotype relative risk \n";
	print "-d: dominance (multiplicative, additive, dominant, recessive) \n";
}



my $i; my $j; my $freq; my $power; my $genotype_relative_risk; my $prev; my $ncases; my $ncontrols; my $alpha;
my $tech_file_nrows; my $s;
my $file_anc_der_raf;
my $line_anc_der_raf; my $line; my @temp=(); my @EAS_DAF=(); my @AMR_DAF=(); my @AFR_DAF=(); my @EUR_DAF=(); my @SAS_DAF=(); my $power_derived; my $power_ancestral;


open (tech_fh_in, "<". $technology);
$tech_file_nrows++ while (<tech_fh_in>);
close (tech_fh_in);


my $R_script;
$R_script = "$working_directory/scripts/weighted_frequency.r";


#********************#
#All parameters fixed#
#********************#

if (defined $cases && defined $controls && defined $significance && defined $prevalence && defined $grr) {
	
	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Reps", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {	
		$genotype_relative_risk = $grr;
		$prev = $prevalence;
		$ncases = $cases;
		$ncontrols = $controls;
		$alpha = $significance;
		
		gas_power();
		
		`Rscript $R_script`;
	
		$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
		open (fh_anc_der_raf, "<". $file_anc_der_raf);
		$line_anc_der_raf = <fh_anc_der_raf>;
		print $fh_output $i+1, "\t", $line_anc_der_raf;
		
		print $i+1, "\n";
	}
	close($fh_output);
}




my $cases_file_nrows; my $controls_file_nrows;
my $significance_file_nrows; my $prevalence_file_nrows; my $grr_file_nrows;

#**********************#
#Number of cases varied#
#**********************#

if (defined $cases_file && defined $controls) {
	my @temp=(); my @cases_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (cases_fh_in, "<". $cases_file);
	$cases_file_nrows++ while (<cases_fh_in>);
	close (cases_fh_in);

	open (cases_fh_in, "<". $cases_file);
	for ($s=0; $s<$cases_file_nrows; $s++) {
		$cases = <cases_fh_in>;
		@temp = split("\n", $cases);
		$cases_list[$s] = $temp[0];
	}
	close (cases_fh_in);

	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Cases_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$cases_file_nrows; $s++) {
			$genotype_relative_risk = $grr;
			$prev = $prevalence;
			$ncases = $cases_list[$s];
			$ncontrols = $controls;
			$alpha = $significance;

			gas_power();
			
			`Rscript $R_script`;
			
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $ncases, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close($fh_output);
}

#*************************#
#Number of controls varied#
#*************************#

if (defined $controls_file && defined $cases) {
	my @temp=(); my @controls_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (controls_fh_in, "<". $controls_file);
	$controls_file_nrows++ while (<controls_fh_in>);
	close (controls_fh_in);
	
	open (controls_fh_in, "<". $controls_file);
	for ($s=0; $s<$controls_file_nrows; $s++) {
		$controls = <controls_fh_in>;
		@temp = split("\n", $controls);
		$controls_list[$s] = $temp[0];
	}
	close (controls_fh_in);
	
	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Controls_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$controls_file_nrows; $s++) {
			$genotype_relative_risk = $grr;
			$prev = $prevalence;
			$ncases = $cases;
			$ncontrols = $controls_list[$s];
			$alpha = $significance;
			
			gas_power();
			
			`Rscript $R_script`;
			
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $ncontrols, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close ($fh_output);
}

#***********************************#
#Number of cases and controls varied#
#***********************************#

if (defined $cases_file && defined $controls_file) {
	my @temp=(); my @cases_list=(); my @controls_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (cases_fh_in, "<". $cases_file);
	$cases_file_nrows++ while (<cases_fh_in>);
	close (cases_fh_in);
	
	open (cases_fh_in, "<". $cases_file);
	for ($s=0; $s<$cases_file_nrows; $s++) {
		$cases = <cases_fh_in>;
		@temp = split("\n", $cases);
		$cases_list[$s] = $temp[0];
	}
	close (cases_fh_in);
	
	open (controls_fh_in, "<". $controls_file);
	$controls_file_nrows++ while (<controls_fh_in>);
	close (controls_fh_in);
	
	open (controls_fh_in, "<". $controls_file);
	for ($s=0; $s<$controls_file_nrows; $s++) {
		$controls = <controls_fh_in>;
		@temp = split("\n", $controls);
		$controls_list[$s] = $temp[0];
	}
	close (controls_fh_in);

	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Cases_sample_size", "\t", "Controls_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$cases_file_nrows; $s++) {
			$genotype_relative_risk = $grr;
			$prev = $prevalence;
			$ncases = $cases_list[$s];
			$ncontrols = $controls_list[$s];
			$alpha = $significance;
			
			gas_power();
			
			`Rscript $R_script`;
	
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $ncases, "\t", $ncontrols, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close ($fh_output);
}
	
#**************#
#P-Value varied#
#**************#

if (defined $significance_file) {
	my @temp=(); my @significance_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (sig_fh_in, "<". $significance_file);
	$significance_file_nrows++ while (<sig_fh_in>);
	close (sig_fh_in);
	
	open (sig_fh_in, "<". $significance_file);
	for ($s=0; $s<$significance_file_nrows; $s++) {
		$significance = <sig_fh_in>;
		@temp = split("\n", $significance);
		$significance_list[$s] = $temp[0];
	}
	close (sig_fh_in);
	
	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Cases_sample_size", "\t", "Controls_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$significance_file_nrows; $s++) {
			$genotype_relative_risk = $grr;
			$prev = $prevalence;
			$ncases = $cases;
			$ncontrols = $controls;
			$alpha = $significance_list[$s];
		
			gas_power();
			
			`Rscript $R_script`;
			
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $alpha, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close ($fh_output);
}

#**********#
#Prevalence#
#**********#

if (defined $prevalence_file) {
	my @temp=(); my @prevalence_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (prev_fh_in, "<". $prevalence_file);
	$prevalence_file_nrows++ while (<prev_fh_in>);
	close (prev_fh_in);
	
	open (prev_fh_in, "<". $prevalence_file);
	for ($s=0; $s<$prevalence_file_nrows; $s++) {
		$prevalence = <prev_fh_in>;
		@temp = split("\n", $prevalence);
		$prevalence_list[$s] = $temp[0];
	}
	close (prev_fh_in);
	
	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Cases_sample_size", "\t", "Controls_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$prevalence_file_nrows; $s++) {
			$genotype_relative_risk = $grr;
			$prev = $prevalence_list[$s];
			$ncases = $cases;
			$ncontrols = $controls;
			$alpha = $significance;
			
			gas_power();
			
			`Rscript $R_script`;
			
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $prev, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close ($fh_output);		
}

#***************#
#GRR (OR) varied#
#***************#

if (defined $grr_file) {
	my @temp=(); my @grr_list=(); my $file_anc_der_raf; my $line_anc_der_raf;
	
	open (grr_fh_in, "<". $grr_file);
	$grr_file_nrows++ while (<grr_fh_in>);
	close (grr_fh_in);
	
	open (grr_fh_in, "<". $grr_file);
	for ($s=0; $s<$grr_file_nrows; $s++) {
		$grr = <grr_fh_in>;
		@temp = split("\n", $grr);
		$grr_list[$s] = $temp[0];
	}
	close (grr_fh_in);
	
	open my $fh_output, '>', "$working_directory/simulation.txt" or die;
	print $fh_output "Cases_sample_size", "\t", "Controls_sample_size", "\t", "AFR_ANC", "\t", "AMR_ANC", "\t", "EAS_ANC", "\t", "EUR_ANC", "\t", "SAS_ANC", "\t", "AFR_DER", "\t", "AMR_DER", "\t", "EAS_DER", "\t", "EUR_DER", "\t", "SAS_DER", "\n";
	
	for ($i=0; $i<$replicates; $i++) {
		for ($s=0; $s<$grr_file_nrows; $s++) {
			$genotype_relative_risk = $grr_list[$s];
			$prev = $prevalence;
			$ncases = $cases;
			$ncontrols = $controls;
			$alpha = $significance;
			
			gas_power();
			
			`Rscript $R_script`;
			
			$file_anc_der_raf = "$working_directory/ANC_DER_RAF.txt";
			open (fh_anc_der_raf, "<". $file_anc_der_raf);
			$line_anc_der_raf = <fh_anc_der_raf>;
			print $fh_output $genotype_relative_risk, "\t", $line_anc_der_raf;
		}
		print $i+1, "\n";
	}
	close ($fh_output);		
}



#************#
#Sub routines#
#************#

#*******************************************************************#
#Calculating power for derived and ancestral risk allele frequencies#
#*******************************************************************#

sub gas_power {
	my $line; my @temp=(); my @EAS_DAF=(); my @AMR_DAF=(); my @AFR_DAF=(); my @EUR_DAF=(); my @SAS_DAF=(); my $j; my $power_derived; my $power_ancestral;
	open my $fh_out_derived, '>', "$working_directory/gas_power_derived.txt" or die;
	print $fh_out_derived "POWER", "\n";
	open (tech_fh_in, "<". "$technology");
	$line = <tech_fh_in>;
	for ($j=1; $j<$tech_file_nrows; $j++) {
		$line = <tech_fh_in>;
		@temp = split("	", $line);
		$EAS_DAF[$i] = $temp[6];
		$AMR_DAF[$i] = $temp[7];
		$AFR_DAF[$i] = $temp[8];
		$EUR_DAF[$i] = $temp[9];
		$SAS_DAF[$i] = $temp[10];

		if ($study_pop eq "eas" || $study_pop eq "EAS" || $study_pop eq "Eas") {
			$freq = $EAS_DAF[$i];
		}
		elsif ($study_pop eq "amr" || $study_pop eq "AMR" || $study_pop eq "Amr") {
			$freq = $AMR_DAF[$i];
		}
		elsif ($study_pop eq "afr" || $study_pop eq "AFR" || $study_pop eq "Afr") {
			$freq = $AFR_DAF[$i];
		}
		elsif ($study_pop eq "eur" || $study_pop eq "EUR" || $study_pop eq "Eur") {
			$freq = $EUR_DAF[$i];
		}
		elsif ($study_pop eq "sas" || $study_pop eq "SAS" || $study_pop eq "Sas") {
			$freq = $SAS_DAF[$i];
		}
			
		if ($freq == 0 || $freq == 1) {
			$power = 0;
			$power_derived = $power;
			print $fh_out_derived $power_derived, "\n";
		}
		else {
			calculate ($ncases, $ncontrols, $freq, $genotype_relative_risk, $prev, $alpha);
			$power_derived = $power - 1;
			print $fh_out_derived $power_derived, "\n";
		}
	}
	close (tech_fh_in);
	close ($fh_out_derived);
		
	open my $fh_out_ancestral, '>', "$working_directory/gas_power_ancestral.txt" or die;
	print $fh_out_ancestral "POWER", "\n";
	open (tech_fh_in, "<". "$technology");
	$line = <tech_fh_in>;
	for ($j=1; $j<$tech_file_nrows; $j++) {
		$line = <tech_fh_in>;
		@temp = split("	", $line);
		$EAS_DAF[$i] = 1 - $temp[6];
		$AMR_DAF[$i] = 1 - $temp[7];
		$AFR_DAF[$i] = 1 - $temp[8];
		$EUR_DAF[$i] = 1 - $temp[9];
		$SAS_DAF[$i] = 1 - $temp[10];
			
		if ($study_pop eq "eas" || $study_pop eq "EAS" || $study_pop eq "Eas") {
			$freq = $EAS_DAF[$i];
		}
		elsif ($study_pop eq "amr" || $study_pop eq "AMR" || $study_pop eq "Amr") {
			$freq = $AMR_DAF[$i];
		}
		elsif ($study_pop eq "afr" || $study_pop eq "AFR" || $study_pop eq "Afr") {
			$freq = $AFR_DAF[$i];
		}
		elsif ($study_pop eq "eur" || $study_pop eq "EUR" || $study_pop eq "Eur") {
			$freq = $EUR_DAF[$i];
		}
		elsif ($study_pop eq "sas" || $study_pop eq "SAS" || $study_pop eq "Sas") {
			$freq = $SAS_DAF[$i];
		}
			
		if ($freq == 0 || $freq == 1) {
			$power = 0;
			$power_ancestral = $power;
			print $fh_out_ancestral $power_ancestral, "\n";
		}
		else {
			calculate ($ncases, $ncontrols, $freq, $genotype_relative_risk, $prev, $alpha);
			$power_ancestral = $power - 1;
			print $fh_out_ancestral $power_ancestral, "\n";
		}
	}
	close (tech_fh_in);
	close ($fh_out_ancestral);

}



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


sub calculate {

	my @pp=(); my $aa_freq; my $ab_freq; my $bb_freq; my @f=();
	my @additive=(); my @multiplicative=(); my @dominant=(); my @recessive=();
	my $scale; my $aa_prob; my $ab_prob; my $bb_prob;
	my $C; my $pcases; my $pcontrols; my $vcases; my $vcontrols; my $ncp;

	
	@pp = ($freq * $freq, 2 * $freq * (1 - $freq), (1. - $freq) * (1. - $freq));

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
	$power = ndist( -1 * $C - $ncp, 'false') + ndist($C - $ncp, 'true');
	return $power;
}