#!/usr/bin/perl 
use strict;
use Data::Dumper;

#################################################################################################################################
#Project: Generating FBA model from source .xml raw data.																							#
#@Author: Madhura Adavkar#
#This perl CODE is basically XML file parser; which is taking .xls derived .csv FILE AS COMMAND LINE ARGUMENT and list out all reaction names ,stoichiometry and genes and#
# generates FBA modeling compatible LIM model as per LIM conventions.																#
##################################################################################################################################

###################################################FUNCTION DEFINITIONS###########################################################
sub Clear{
	my $catch = $_[0];
	$catch=~s/^\s+//g;
	$catch=~s/\s+$//g;
	return($catch);
}

#This is function for getting all unique metabolites
sub Get_uniq_met{
	my $genearr_ref = shift;
	my @genearr = @{$genearr_ref};
	my %Unique_genes;
	my @uniq_genearr=();

	foreach my $gg(@genearr){
		$gg =Clear($gg);	
		if(exists $Unique_genes{$gg}){next;}
		else{
			$Unique_genes{$gg}=1;
			push(@uniq_genearr,$gg);
		}			
	}
	
	my $before = scalar(@genearr);
	my $after = scalar(@uniq_genearr);
	
	#print "Previous number of genes: $before\t After removing repeats, unique number of genes: $after \n";
	return(\@uniq_genearr);
}

#This function will list out all metabolites in reactions. Giving all unique metabolites in sorted order.
sub Eval_Metabolites{
	my $Ref_got = shift;
	my @Sto_arr = @{$Ref_got};
	my @Compo_arr =();
	
	foreach my $Stt(@Sto_arr){
		$Stt = Clear($Stt);
		my $sep ;
		
		if($Stt =~/\<\-\>/){
			$sep = "\<\-\>";
		}elsif($Stt =~/\<\-/){
			$sep = "\<\-";					
		}elsif($Stt =~/\-\>/){
			$sep = "\-\>";				
		}else{
		}
		
		my @sp_stt = split($sep,$Stt);
		foreach my $pseudo_met(@sp_stt){
			$pseudo_met = Clear($pseudo_met);
			if($pseudo_met =~/\+/){
				my @sp2 = split(/\+/,$pseudo_met);
				foreach my $meta(@sp2){
					$meta = Clear($meta);
					if($meta =~/((\d+)\.?(\d+)?)\*/){
							$meta =~s/((\d+)\.?(\d+)?)\*//g;
							$meta = Clear($meta);
							push(@Compo_arr , $meta);				
					}else{
						push(@Compo_arr , $meta);
					}
					
				}		
			}else{
			#	push(@Compo_arr , $pseudo_met);
				if($pseudo_met =~/((\d+)\.?(\d+)?)\*/){
						$pseudo_met =~s/((\d+)\.?(\d+)?)\*//g;
						$pseudo_met= Clear($pseudo_met);
						#print "$meta for $pseudo_met \n";
						push(@Compo_arr , $pseudo_met);				
				}else{
					push(@Compo_arr , $pseudo_met);
				}				
			}
		}
					
	}
	
	#Now taking unique components list	
	my $Final_ref = Get_uniq_met(\@Compo_arr);
	my @Final_arr = @{$Final_ref};
	
	my @Uniq_Compo_arr = sort(@Final_arr);			
	return(\@Uniq_Compo_arr);
}

#This is function for evaluating COMPONENTS discription provided in .csv from .xls sheet
sub Get_Component_discp{
	#my $csv_file2 = "msb4100155-s3_Ecoli_Components_discription.csv";
	my $csv_file2 = $ARGV[1];
	
	print "CSV File refereed for Components discription : $csv_file2\n";
	print "****************************************************\n";
	
	my %HASH_Comp_Disc = {};
	
	#Reading csv file(! deliminator) and parsing it
	open(CSV2,$csv_file2)|| die "$!";
	while(<CSV2>){
		my $line = $_;
		chomp($line);	
		#print "$line \n";
		$line = Clear($line);
		if($line eq ""){
			next;
		}
		if($line =~/^Abbreviation/i){
			next;
		}
		my @arr_line2 = split(/\!/,$line);
		my $Met_name = Clear($arr_line2[0]);
		my $Disp = Clear($arr_line2[-1]);
		
		
		if($Met_name=~/\b\-\b/){
			my $Bef = $Met_name;
			$Met_name =~s/\b\-\b/\_/g;
			$Met_name =~s/\[/\_/g;
			$Met_name =~s/\]/ /g;
		}
		
		$Met_name = uc($Met_name);
		
		#print "$Met_name ******* $Disp \n";
		$HASH_Comp_Disc{$Met_name} = $Disp;		
	}		
	close(CSV2);		
	return(\%HASH_Comp_Disc);
}
############################################################################################################

##############################################################################################################
#Taking inputs: 1st command line argument is .csv for LIM model and 2nd arg is components discription file
#Taking command line argument

#my $csv_file = "msb4100155-s3_Ecoli_Mg1655Model_RXN_Listing_Final.csv";
#my $csv_file2 = "msb4100155-s3_Ecoli_Components_discription.csv";

my $csv_file = $ARGV[0];
print "********************INPUTS********************************\n";
print "CSV File referred for LIM model reactions: $csv_file\n";

my $Num_reactions = 0;

#Reading csv file(! deliminator) and parsing it ;ine by line
open(CSV,$csv_file)|| die "$!";

#Definying output file filehandle
my $Output = "Final_Output_EcoliMG1655_LIM_model.lim";		#This is output .lim file (FBA model with LIM model conventions.)
open(OUTLIM ,">$Output")|| die "$!";
open(TRANSPORTRXN,">Final_Transport_rxns_listing_Ecoli.txt")|| die "$!";		#This output file list out all transport reactions
open(RXNGENEASSOC,">Final_Reactions_Gene_Association_EcoliMG1655.txt")|| die "$!";  		##This output file list out all genes associated with respective reactions
open(RXNCSV,">Final_CSV_Format_Ecoli.csv")|| die "$!";				##This output file list out .csv file
open(COMPONENTSCSV,">Final_Components_Ecoli.csv")|| die "$!";		#This output file list out all metabolites
open(EXTERNALSCSV,">Final_Externals_Ecoli.csv")|| die "$!";			##This output file list out all externals

############Inits
my %HASH_RXN_Stoich = {};
my @Arr_Stoich = ();
my @Arr_Subsystem = ();

print RXNCSV "Reaction name\t!\tStoichiometry\t!\tEnzyme discription\t!\tCompartment\t!\tSubsystem(Pathway)\t!\tGene\n";
print OUTLIM "## REACTIONS\n";
while(<CSV>){
	my $line = $_;
	chomp($line);	
	#print "$line \n";
	$line = Clear($line);
	if($line eq ""){
			next;
	}
	if($line =~/^Abbreviation/i){
		next;
	}
	
	my @arr_line = split(/\!\!?/,$line);

	my $Rxn_name = Clear($arr_line[0]);
	my $Stoich = Clear($arr_line[1]);
	my $Enzyme_disc_Name = Clear($arr_line[2]);
	my $Subsystem = Clear($arr_line[3]);
	my $Gene = Clear($arr_line[-1]);
	
	#my $Protein = Clear($arr_line[-2]);
	#my $Subsystem = Clear($arr_line[-3]);
	#my $Protein = Clear($arr_line[-2]);
	#my $Gene = Clear($arr_line[-1]);
	my $Compartment ;
			
	
	#Replacing hypen by underscore and compartment name as _c or _h or _e
	if($Rxn_name =~/\b\-\b/){
		$Rxn_name =~s/\b\-\b/\_/g;
	}
	#(e)
	if($Rxn_name =~/\(e\)/){
		my $Bef = $Rxn_name;
		$Rxn_name =~s/\(e\)/\_e/g;
	}	
	
	if($Stoich =~/^\[C\]\s*\:/i){	
		$Compartment = "Cytosol";
		$Stoich =~s/\[C\]\s*\://i;		
	}elsif($Stoich =~/^\[P\]\s*\:/i){
		$Compartment = "Periplasmic space";
		$Stoich =~s/\[P\]\s*\://i;			
	}elsif($Stoich =~/^\[E\]\s*\:/i){
		$Compartment = "Extracellular space";
		$Stoich =~s/\[E\]\s*\://i;	
	}else{
		$Compartment = "Intercompartmental";			
	}
	
	
	#Making Stoichmetry uppercase
	$Stoich = uc($Stoich);	
	
	if($Stoich =~/\[/){
		$Stoich =~s/\[/\_/g;		
		$Stoich =~s/\]/ /g;
	}	
	
	#Making Stoichmetry in terms of number of moles with star between digit and component name
	if($Stoich =~/(\((\d+)\.?(\d+)?\))/){
		$Stoich =~s /\(//g;
		$Stoich =~s /\)/\*/g;
	}	
	
	#Making reaction direction 
	if($Stoich =~/\-\-\>/){			
		$Stoich =~s/\-\-\>/->/;		
	}elsif($Stoich =~/\<\-\-/){
		$Stoich =~s/\<\-\-/<-/;			
	}elsif($Stoich =~/\<\=\=\>/){
		$Stoich =~s/\<\=\=\>/<->/;		
	}else{
	}
	
	#Replacing hypen by underscore and compartment name as _c or _h or _e
	if($Stoich =~/\b\-\b/){
		my $Bef = $Stoich;
		$Stoich =~s/\b\-\b/\_/g;
	}
	
	$Num_reactions++ ;
	
	
	my $Bef = $Stoich;
	$Bef = Clear($Bef);
	#For taking out compartment name from component name
	if($Stoich =~/\[C\]/gi){		
		$Stoich =~s /\[C\]/\_C/gi;		
	}
	if($Stoich =~/\[E\]/gi){		
		$Stoich =~s /\[E\]/\_E/gi;		
	}
	
	$Stoich = Clear($Stoich);
	
	if($Gene =~/(b|s)\d\d\d\d/){		
		if($Gene =~/and/){
			$Gene =~s/and/&&/g;			
		} 
		if($Gene =~/or/){
			$Gene =~s/or/||/g;		
		} 
		print RXNGENEASSOC "$Rxn_name\t:\t$Gene\n";	
		
	}elsif($Gene =~/(\d+)\.(\d+)\.(\d+).(\d+)/){
		print RXNGENEASSOC "$Rxn_name\t:\t$Gene\n";
	
	}else{
		$Gene = "Unknown";
		print RXNGENEASSOC "$Rxn_name\t:\t$Gene\n";
	}	
	

	push(@Arr_Stoich,$Stoich);
	
	print RXNCSV "$Rxn_name\t!\t$Stoich\t!\t$Enzyme_disc_Name\t!\t$Compartment\t!\t$Subsystem\t!\t$Gene\n";
	print OUTLIM "$Rxn_name\t:\t$Stoich\t!\t$Enzyme_disc_Name\t!\t$Compartment\t!\t$Subsystem\t!\t$Gene\n";
	
	$HASH_RXN_Stoich{$Rxn_name} = $Stoich ;
	
	#Evaluating Transporter reactions
	if($Rxn_name =~/tex(i)?$/){
		print TRANSPORTRXN 	"$Rxn_name\t:\t$Stoich\t!\t$Enzyme_disc_Name\t!\t$Compartment\t!\t$Subsystem\t!\t$Gene\n";
	}elsif($Rxn_name =~/ex(i)?$/){
		print TRANSPORTRXN 	"$Rxn_name\t:\t$Stoich\t!\t$Enzyme_disc_Name\t!\t$Compartment\t!\t$Subsystem\t!\t$Gene\n";
	}elsif($Rxn_name =~/pp(i)?$/){
		print TRANSPORTRXN 	"$Rxn_name\t:\t$Stoich\t!\t$Enzyme_disc_Name\t!\t$Compartment\t!\t$Subsystem\t!\t$Gene\n";
	}
	
	push(@Arr_Subsystem,$Subsystem);
}

print OUTLIM "## END REACTIONS\n";

#For Objective Function section
print OUTLIM "## MAXIMISE\n";

print OUTLIM "## END MAXIMISE\n";

#For Constraints section
print OUTLIM "## INEQUALITIES\n";
foreach my $Rx(keys(%HASH_RXN_Stoich)){
	$Rx = Clear($Rx);
	if($Rx =~/^HASH\(/){
		next;
	}

	my $Sto = $HASH_RXN_Stoich{$Rx};	
	my $Constraints ;
	
	#Making reaction direction 
	if($Sto =~/\<\-\>/){			
		$Constraints = "[\t-999999\t,\t999999\t]" ;
	}elsif($Sto =~/\<\-/){
		$Constraints = "[-999999\t,\t0\t]" ;		
	}elsif($Sto =~/\-\>/){
		$Constraints = "[\t0\t,\t999999\t]" ;
	}else{
	}
	print OUTLIM "$Rx\t\=\t$Constraints\t\!\t$Sto\n" ;
		
}
print OUTLIM "## END INEQUALITIES\n";
print OUTLIM "\n";

##################################################################################
my $Ref_Results_list = &Eval_Metabolites(\@Arr_Stoich);
my @Final_met_arr = @{$Ref_Results_list};

#For Components_section
print OUTLIM "## COMPONENTS\n";
my $ref_hash = &Get_Component_discp();
my %Comp_Dis_Got = %{$ref_hash};
#print Dumper\%Comp_Dis_Got;

my $c = 0;
my @Externals_arr =();
foreach my $Comp(@Final_met_arr){
	$Comp = Clear($Comp);
	#Checking if component disciption exits or not
	if(exists($Comp_Dis_Got{$Comp})){
		#print "Exists disp for $Comp as $Comp_Dis_Got{$Comp} \n";
		$c++;
		if($Comp =~/_E$/){
			push(@Externals_arr,$Comp);
		}else{
			print OUTLIM "$Comp ! $Comp_Dis_Got{$Comp}\n";	
			print COMPONENTSCSV "$Comp ! $Comp_Dis_Got{$Comp}\n";
		}		
		
	}else{
		if($Comp =~/_E$/){			
			push(@Externals_arr,$Comp);
		}else{
			print OUTLIM "$Comp ! $Comp_Dis_Got{$Comp}\n";	
			print COMPONENTSCSV "$Comp ! $Comp_Dis_Got{$Comp}\n";
		}		
	}
}
print OUTLIM "## END COMPONENTS\n";

print OUTLIM "\n";

print OUTLIM "## EXTERNALS\n";
foreach my $ext(@Externals_arr){
	$ext = Clear($ext);
	if(exists($Comp_Dis_Got{$ext})){
		print OUTLIM "$ext ! $Comp_Dis_Got{$ext}\n";	
		print EXTERNALSCSV "$ext ! $Comp_Dis_Got{$ext}\n";
	}else{
		print OUTLIM "$ext ! \n";
		print EXTERNALSCSV "$ext ! \n";
	}
}
print OUTLIM "## END EXTERNALS\n";

print "\n";
######################################################################################
print "**********************OUTPUT******************************\n";
print "\033[31m Number of evaluated component :". scalar(@Final_met_arr)."\033[0m \n";

close(COMPONENTSCSV);
close(EXTERNALSCSV);
close(RXNCSV);
close(RXNGENEASSOC);
close(TRANSPORTRXN);
close(OUTLIM);
close(CSV);


print "\033[31m Number of reactions in model: $Num_reactions \033[0m \n";

print "*******************OUTPUT FILES ********************************\n";
print "\033[31m 1) E coli Model: Final_Output_EcoliMG1655_LIM_model.lim \033[0m \n";
print "\033[31m 2) Text file listing all transport reactions: Final_Transport_rxns_listing_Ecoli.txt \033[0m \n";
print "\033[31m 3) Text file showing Reactions and genes association:Final_Reactions_Gene_Association_EcoliMG1655.txt \033[0m\n";
print "\033[31m 4) LIM model reactions in csv format: Final_CSV_Format_Ecoli.csv \033[0m \n";
print "\033[31m 5) E coli Components List:Final_Components_Ecoli.csv \033[0m\n";
print "\033[31m 6) E coli Externals List: Final_Externals_Ecoli.csv \033[0m\n";
print "\033[31m 7) E coli Pathway List(Subsystems): Pathways_Listing.txt \033[0m\n";
print "****************************************************\n";
#Now taking unique Pathways list	

my $Final_ref_sub = Get_uniq_met(\@Arr_Subsystem);
my @Final_arr_sub = @{$Final_ref_sub};

open(PATHWAYS,">Pathways_listing.txt")|| die "$!";
my @Uniq_Pathway_arr = sort(@Final_arr_sub);		
foreach my $path(@Uniq_Pathway_arr){
	$path = Clear($path);
	#print "$path\n";
	print PATHWAYS "$path\n";
}
close(PATHWAYS);
#print Dumper\%HASH_RXN_Stoich;

