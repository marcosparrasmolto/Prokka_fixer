#!/usr/bin/perl

use strict;
#use warnings;

#############################################################################################################################################################################
###                                                                                                                                                                       ###
###		Prokka_parser. Get custom database information into Prokka gff filefield																						  ###
###                                                                                                                                                                       ###
###                                                                                                                                                                       ###
#############################################################################################################################################################################

my $dir_blast="./";

my %list_orf_blast_pos=undef;
my %list_orf_blast=undef;
my %list_comb_entrada=undef;

open (entrada_blast,$dir_blast."aux_file_prodigal_global_headers.txt")  || die "File could not be opened\n";

while(<entrada_blast>) ##Information abour each ORF sequence are saved.
{
	chomp $_;
	
	my @line_blast=split(/ \# /,$_);
	$list_orf_blast_pos{substr($line_blast[0], 1)}=$line_blast[1].";".$line_blast[2];
	#print substr($line_blast[0], 1)." ".$list_orf_blast_pos{substr($line_blast[0], 1)}."\n";
	#print $line_blast[0]."_".($list_orf_blast{$line_blast[0]})."\n";	
}



close entrada_blast;

open (entrada_blast,$dir_blast."All_results_hit.txt") || die "File could not be opened\n"; ##Best blast hit for each ORF file are read and information is saved

while(<entrada_blast>)
{
	chomp $_;
	
	my @line_blast=split(/\t/,$_);
	$list_orf_blast{$line_blast[0]}=$line_blast[1].";".$line_blast[13];
	#print $line_blast[0]." ".($list_orf_blast{$line_blast[0]})."\n";	
	
	
}


foreach my $key (keys %list_orf_blast) ##Information about blast information according to ORF information, such as starting and ending position is saved into arrays
{
		# do whatever you want with $key and $value here ...
		#my $value = $list_orf_blast_pos{$key};
		#print "  $key costs $value\n";
		
		my $i=rindex($key,"_");
		my $a=substr($key,0,$i);
		#print $key." ".$a."\n";
		
		
		$list_comb_entrada{$a.";".$list_orf_blast_pos{$key}}=$list_orf_blast{$key};
		#print $key." ".$a.";".$list_orf_blast_pos{$key}." ".$list_orf_blast{$key}."\n";
}	





#=pod
my $dir="./Output";

my @files = glob( $dir . '/*' );

my %list_orf=undef;
my $flag_gene=0;
my $flag_name=0;
my $flag_product=0;
my $flag_note=0;
my $back_product;

for(@files) ##For each gff file information about ORFs with new blast information from custom database we open and edit them according to the new information
{
	my @files2 = glob( $_ . '/*' );
	my $dir_guardado=$_;
	
	
	
	for(@files2)
	{		
		if($_ =~ /gff/)
		{
			$dir_guardado=(split(/.contigs.fa/,$dir_guardado))[0];
			
			open (entrada,$_) || die "File could not be opened\n";
			open (salida,">".$dir_guardado."_new_prokka.gff");

			while(<entrada>)
			{
				chomp $_;
				if($_ =~ /ID=/) ##Each time we find an ID field we check if we have new blast information or not
				{			
					my @line=split(/\t/,$_);
					
					if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]})
					{
						
						my @line_ID=split(/\;/,$line[8]);
						$flag_gene=0;
						$flag_product=0;
						$flag_name=0;
						$flag_note=0;	
						$line[8]=undef;
						$back_product=undef;
						
						foreach my $linea_ID_ver (@line_ID) ##We will check each field: gene, name, product and note, and will replace the previous information for the new one in the case we got a blast result from the custom database
						{
							if($linea_ID_ver =~ /gene=/)
							{
								#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
								#{
								#	$linea_ID_ver="gene=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[1];
								
								#}else
								#{
									$linea_ID_ver="gene=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[1];
								#}
									
								$flag_gene=1;
							}if($linea_ID_ver =~ /Name=/)
							{
								#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
								#{
								#	$linea_ID_ver="Name=".(split(/\_/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[2];
								
								#}else
								#{
									$linea_ID_ver="Name=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[1];
								#}
									
								$flag_name=1;
							}elsif($linea_ID_ver =~ /product=/)
							{
								$back_product=(split(/product=/,$linea_ID_ver))[1];


								#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
								#{
								#	$linea_ID_ver="product=";
								#	my @product_line_=split(/_/,(split(/\;/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[0]);
								#	@product_line_=@product_line_[3..$#product_line_];
								#	$linea_ID_ver.="$_ " for @product_line_;
								#	chop($linea_ID_ver);
								#}else
								#{
									$linea_ID_ver="product=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[2];

									#$linea_ID_ver="product=";
									#my @product_line_=split(/?/,(split(/\;/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[0]);
									#@product_line_=@product_line_[2..$#product_line_];
									#$linea_ID_ver.="$_ " for @product_line_;
									#chop($linea_ID_ver);
								#}

								$flag_product=1;
							}elsif($linea_ID_ver =~ /note=/)
							{
								$linea_ID_ver="note=".$back_product."-Annotation_modified";
								$flag_note=1;
							}
							
							$line[8].=$linea_ID_ver.";";
						}
						
						
						
						if($flag_gene==0) ##If the original field doesnt contain a field, we will create the entry for it
						{
							#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
							#{
							#	$line[8].="gene=".(split(/\_/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[2].";";
							#}else
							#{
								$line[8].="gene=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[1].";";
							#}
						}
						if($flag_name==0)
						{
							#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
							#{
							#	$line[8].="Name=".(split(/\_/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[2].";";
							#}else
							#{
								$line[8].="Name=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[1].";";
							#}
						}

						if($flag_product==0)
						{
								#if($list_comb_entrada{$line[0].";".$line[3].";".$line[4]} =~ /NZ_/)
								#{
								#	$line[8].="product=";
								#	my @product_line_=split(/_/,(split(/\;/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[0]);
								#	@product_line_=@product_line_[3..$#product_line_];
								#	$line[8].="$_ " for @product_line_;
								#	chop($line[8]);
								
								#}else
								#{
									$line[8].="product=".(split(/\?/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[2];

									#$line[8].="product=";
									#my @product_line_=split(/_/,(split(/\;/,$list_comb_entrada{$line[0].";".$line[3].";".$line[4]}))[0]);
									#@product_line_=@product_line_[2..$#product_line_];
									#$line[8].="$_ " for @product_line_;
									#chop($line[8]);
								#}
						}
						if($flag_note==0)
						{
							$line[8].="note=".$back_product."-Annotation_modified;";
						}
						
						chop($line[8]);

						print salida "$_\t" for @line; ##Information is saved line by line
						print salida "\n";
					
					}else
					{
						print salida "$_\t" for @line;
						print salida "\n";
					}					
				}else
				{	
					print salida $_."\n";
				}
			}
			
			close entrada;
			close salida;
		}
	}
	
	%list_orf=undef;
}	
