#!/usr/bin/perl -w
#Authors: Alina Saiakhova (ars51@case.edu) and Olivia Corradin (ogc@case.edu)
BEGIN {$^W = 0}
use Cwd;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;
use List::Util qw[min max];
use File::Copy;

=head1 NAME

    sample - Using PreSTIGE_pipeline.pl

    
=head1 SYNOPSIS

Example: ./PreSTIGE_pipeline.pl  -p /full/path/to/cell_line_H3K4me1_peak_file.txt -w /full/path/to/wig  -n cell_line_name -e /full/path/to/cell_line_FPKM_file.txt

    Options:
    
    -h|help           Prints this message
    
    -peaks|p                Full path to new cell line peak file
    
    -wig|w                Full path to all-chromosome wig track for new cell line
    
    -name|n                Cell line name
    -expr|e                Full path to RNA-seq or microarray expression data for new cell line. The file must be in this format <gene name><TAB or space><expression>.
    
    -m                (Optional) Set -m to specify microarray expression data. Otherwise, RNA-seq data is assumed to be used.
    
    -gl               (Optional) Merge by gene locus (from cufflinks output)
    -dbkey            (Optional) Genome assembly (default is hg18)
    -acpeaks           Full path to H3K27ac peaks file
   -hs                Get high stringency predictions
   -t,-thread                 Number of cores to use with Heatmap Script (maximum is 4)



=cut


my %options=('dbkey'=>'hg18');
 if (exists($options{'help'}) or @ARGV<3){pod2usage(1);}
GetOptions (\%options, 'help|h|?' , 'peaks|p=s','wig|w=s','name|n=s','expr|e=s','m','gl', 'dbkey=s','acpeaks=s','hs', 'thread|t=i');

my $working_dir=getcwd;
my $log_file="log.txt"; open(LOG, ">$log_file") || print "Can't open log file\n"; ##log file code

##validate COMMAND-LINE PARAMETERS
errorHandler("","Input parameter error:Peak file parameter not specified") if !defined($options{'peaks'});
errorHandler("","Input parameter error:Wig file parameter not specified") if !defined($options{'wig'});
errorHandler("","Input parameter error:Sample name parameter not specified") if !defined($options{'name'});
errorHandler("","Input parameter error:User FPKM file parameter not specified") if !defined($options{'expr'});

#INITIALIZE ALL COMMAND-LINE VARIABLES
my $peak_file=$options{'peaks'}; ##new_sample_peaks_file
my $wig_file=$options{'wig'}; ##/home1/ars51/OLIVIA/PRESTIGE_PIPELINE/new_sample_all_chromosome_wig.wig
my $sample_name=$options{'name'}; ##new_sample_name
my $user_fpkm_file=$options{'expr'}; ##expression file
my $genome=$options{dbkey}; ##genome assembly
my $ac_peaks_file=$options{acpeaks} if defined($options{acpeaks});
my $num_cores=1;
$num_cores=$options{thread} if defined($options{thread});



##INTERNAL VARIABLES (SUBJECT TO CHANGE)
my $scripts="full_path_to_scripts"; ##this WILL change between evolution and mendel
my $txt_dir="full_path_to_txt_files";  ##this WILL change between evolution and mendel
my $db_wig_dir="full_path_to_wigs/$genome/"; ##this WILL change between evolution and mendel
my $db_enh=$txt_dir."/$genome/MANU_enhancer_peaks_all"; ##merged peaks file for all db samples
my $db_ctcf=$txt_dir."/$genome/g40_chr_midpoint.txt"; ##ctcf sites
my $db_d3=$txt_dir."/$genome/Chromosome_start_stop_getcounts.txt"; ##file used to getCounts using distance model
my $tss_for_getcounts=$txt_dir."/$genome/rFLAT_tss_for_getcounts.txt"; ##TSS used for getCounts
my $db_fpkm=$txt_dir."/$genome/FPKM_12lines";  ##FPKM matrix of all db samples
my $c_names=$txt_dir."/$genome/col1_for_getcounts_geneid.txt"; ## tss gene names for getCounts
my $hm_dir="heatmap_script"; ##location where all HM output files will be stored
my $hm_subdir="Heatmap-MinMedianMax_3.1"; ##version of HM script used
my $samples_comma_delim_string="GM12878,H1ES,HepG2,HMEC,HSMM,HUVEC,K562,NHEK,NHLF,HeLa,MCF7,NPC"; ## these are the database samples that will be used
my $chr_for_hm_script="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X";
if($genome eq 'mm9'){$samples_comma_delim_string="cerebellum,embryonic-heart,embryonic-limb,intestine,kidney,liver,MEF,olfactory-bulb,placenta,spleen,testis,thymus"; $chr_for_hm_script="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X"}
##END OF INTERNAL VARIABLES

my @samples_arr=();
if($samples_comma_delim_string =~m/,/){@samples_arr=split(/,/, $samples_comma_delim_string);}else{ push @samples_arr, $samples_comma_delim_string;}


foreach(@samples_arr){
   
   if($sample_name eq $_){
      
      errorHandler("", "$sample_name already exists in the PreSTIGE database and cannot be used. You must provide a different sample name.\n");
      
   }
}


print LOG "Parameters:\npeak_file:$peak_file\nwig_dir:$wig_file\nsample_name:$sample_name\nuser_fpkm_file:$user_fpkm_file\n";

if(! exists($options{'m'})){print LOG "type_of_expr:RNA-Seq expression\n";}else{print LOG "type_of_expr:Microarray data\n";}

###MAIN PIPELINE STARTED
@timeData = localtime(time);
print LOG "Main pipeline started on $genome: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";
 print LOG "Checking user FPKM file format...";
       $user_fpkm_file=checkUserFPKMFileFormat($user_fpkm_file);
print LOG "Done\n";
checkUserPeakFileFormat($peak_file);

##Determine which chromosomes to run script on
print LOG "Determining chromosomes to analyze...\n";
($chr_for_hm_script, $db_enh, $peak_file)=checkIfAllChromsInUserFiles($peak_file, $wig_file, $db_enh);
#($chr_for_hm_script, $db_enh, $peak_file)=($chr_for_hm_script, "all_merged_enhancers.txt", "new_sample_peak_list.txt");
if($chr_for_hm_script == ""){print LOG "No valid chromosomes found in user files. Exiting.\n"; die;}
print LOG "Done\n";



runStep1($peak_file);
runStep2($wig_file, $peak_file, $sample_name);
runSteps3and4($user_fpkm_file,$db_fpkm,$sample_name);
@timeData = localtime(time);
print LOG "Main pipeline finished: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";



#SHARED SUBROUTINES
sub errorHandler{
my $standard_error_string=$_[0];
my $custom_error_string=$_[1];
print LOG "$custom_error_string $standard_error_string\n";
die;
	
}
sub log2 {
	my $n = shift;
	return (log($n)/ log(2));
}
sub shannon_entropy_scoring{
print LOG "Shannon entropy scoring $_[0]...";
    my $qnormed_file=$_[0];
    my $output_file=$_[1];
    open(OUTPUT, ">$output_file");
    open(QNORMED, "<$qnormed_file") or print LOG "Unable to open q-normed file $!\n";   
    while (my $line=<QNORMED>){
	chomp $line;
        if($line=~m/H3K4me1|gene_name/){
         my @liner=split(/\t|\s+/, $line);
        
         print OUTPUT $line."\t";
         
         for(my $i=0; $i<$#liner;$i++){
            if ($liner[$i]=~m/^(.+)\_(H3K4me1|expr)$/){
               print OUTPUT $1."_SE\t";
            }
            
         }
         if ($liner[$#liner]=~m/^(.+)\_(H3K4me1|expr)$/){
         print OUTPUT $1."_SE\n";
         }
        }
        else{
	my @liner=split(/\t|\s+/, $line);
	my $sum=0;
	my @s_arr=();
	
	my @data_only=@liner;
	splice(@data_only, 0, 3);
	
	foreach(@data_only){$sum+=$_;}
	
	@s_arr=map{$_/$sum} @data_only;
	@H_enh_MCF_arr=map{-$_*log2($_)} @s_arr;
	
	my $H_enh_MCF=0;
	foreach(@H_enh_MCF_arr){$H_enh_MCF+=$_;}
	
	my @q_scores=map{$H_enh_MCF-(log2($_))} @s_arr;
	my $last_q_sc=pop @q_scores;
        print OUTPUT (shift @liner)."\t";
        print OUTPUT (shift @liner)."\t";
        print OUTPUT (shift @liner)."\t";
	print OUTPUT map { sprintf("%.4f",$_)."\t" } @liner;
	print OUTPUT map{sprintf("%.4f",$_)."\t"} @q_scores;
	print OUTPUT sprintf("%.4f",$last_q_sc);
	print OUTPUT "\n";
        }
    }

    
    print LOG "Done\n";
}




##STEP 1 SUBROUTINES###
sub runStep1{
print  LOG "Running Step 1\n";
@timeData = localtime(time);
print LOG "Step 1 started: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";

    $peak_file=$_[0];
   # INITIAL ENHANCER-GENE PAIRS IN CTCF DOMAINS

#1. cat all enhancer peaks, merge them, and get their midpoints

errorHandler("","Peak file $peak_file doesn't exist") if (! -f $peak_file);

print LOG "Getting enhancer midpoints...";
system("cat $db_enh $peak_file | sortBed -i stdin|mergeBed -i stdin  | awk '{printf  \"%s\t%d\\n\", \$1, (\$2+\$3)/2}' > all_merged_enhancers_midpoints.txt")==0 or errorHandler($!,"Catting peaks 1 failed");
print LOG "Done\n";
open(FILE, "< all_merged_enhancers_midpoints.txt"); $merged_enh_count++ while (<FILE>);
errorHandler("","User peaks file merged with PreSTIGE merged peaks file must be less than 400,000 lines. Found a total of $merged_enh_count peaks.Exiting pipeline...\n") if $merged_enh_count>=400000;

#2. CTCF approach:
    #run getCountsNewest on: g40(ctcf), rFLAT(tss), enh(enhancer midpoints)
    #append column names
print LOG ("Running getCountsNewest using CTCF-approach...");
system("perl $scripts/getCountsNewest.pl $db_ctcf $tss_for_getcounts all_merged_enhancers_midpoints.txt  ctcf_getcounts_output_tmp.txt 10000000")==0 or errorHandler($!,"getCounts on CTCF failed");
system("paste $c_names ctcf_getcounts_output_tmp.txt > ctcf_getcounts_output.txt");
unlink("ctcf_getcounts_output_tmp.txt");
print LOG "Done\n";
    #parse output to get enhancer-gene pairs in CTCF domains
    print LOG ("Parsing getCountNewest results from CTCF-approach...");
   getEgpFromCountsOutput("ctcf_getcounts_output.txt", "parsed_g40.txt", "g40");
print LOG "Done\n";
#3. DISTANCE approach:
   #run getCountsNewest on: d3(chrom_start_stop file), rFLAT(tss), enh(enhancer midpoints)
    #append column names
    print LOG "Running getCountsNewest using distance-approach...";
    system("perl $scripts/getCountsNewest.pl $db_d3 $tss_for_getcounts all_merged_enhancers_midpoints.txt  d3_getcounts_output_tmp.txt 100000")==0 or errorHandler($!,"getCounts on d3 failed");
system("paste $c_names d3_getcounts_output_tmp.txt > d3_getcounts_output.txt");
unlink("d3_getcounts_output_tmp.txt");
print LOG "Done\n";
    #parse output to get enhancer-gene pairs within 100KB of genes
    print LOG "Parsing getCountNewest results from distance-approach...";
  getEgpFromCountsOutput("d3_getcounts_output.txt", "parsed_d3.txt", "d3");
print LOG "Done\n";
#4. cat enhancer-gene pairs from 2 approaches and label pairs as belonging to CTCF-only, d3-only or both
    print LOG "Catting g40 and d3 files...";
   system("cat parsed_g40.txt parsed_d3.txt > d3_CTCF_catted.txt " )==0 or errorHandler ($!,"Catting d3 and g40 output failed");
print LOG "Done\n";
   print LOG "Getting unique and duplicate gene-enhancer pairs...";
   getDuplicatesBasedOnColumn("d3_CTCF_catted.txt","un_d3_CTCF.txt", "dup_d3_CTCF.txt",1) ;
print LOG "Done\n";
   print LOG "Redistributing gene-enhancer pairs across multiple files...";
sepSharedUniquePairs("un_d3_CTCF.txt", "dup_d3_CTCF.txt", "un_CTCF_only.txt", "un_d3_only.txt", "dup_lab_d3_CTCF.txt" );
print LOG "Done\n";
print  LOG "Get final gene-ehnancer pairs file...";
   system("cat un_CTCF_only.txt un_d3_only.txt dup_lab_d3_CTCF.txt > d3_CTCF_all_lab.txt")==0 or errorHandler ($!,"Catting files failed");
print LOG "Done\n";
    @timeData = localtime(time);
print LOG "Step 1 finished: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";

}
####HELPER SUBROUTINES####
sub sepSharedUniquePairs{
    my $unique_input_file=$_[0];
    my $dup_input_file=$_[1];
    my $output_unique_ctcf_file=$_[2];
    my $output_unique_d3_file=$_[3];
    my $output_shared_file=$_[4];
    
    ##unique
open(UNIQUE_INPUT, "<$unique_input_file") or  die("$!");
open(UNIQUE_CTCF_OUTPUT, ">$output_unique_ctcf_file") or die("$!");
open(UNIQUE_D3_OUTPUT, ">$output_unique_d3_file") or die ("$!");
while(my $line=<UNIQUE_INPUT>){
    chomp $line;
    my @liner=split(/\t|\s+/, $line);
    if($liner[1]=~m/d3/){
        print UNIQUE_D3_OUTPUT "$liner[0]\tNA\td3\n";
    }
    elsif($liner[1]=~m/g40/){
         print UNIQUE_CTCF_OUTPUT "$liner[0]\tg40\tNA\n";
        
    }
    
}
close  UNIQUE_INPUT;
close UNIQUE_D3_OUTPUT;
close UNIQUE_CTCF_OUTPUT;

##shared
open(DUP_INPUT, "<$dup_input_file") or die("$!");
open(SHARED_OUTPUT, ">$output_shared_file") or die ("$!");
my %eh=();
while(my $line=<DUP_INPUT>){
    chomp $line;
    my @liner=split(/\t|\s+/, $line);
    if (! $eh{$liner[0]}++){
    print SHARED_OUTPUT "$liner[0]\tg40\td3\n";
    
    }
    
}

close DUP_INPUT;
close SHARED_OUTPUT;

}



sub getDuplicatesBasedOnColumn{
    my $input_file_name=$_[0];
    my $unique_output_file_name=$_[1];
    my $dup_output_file_name=$_[2];
    my $column=$_[3];
    my %lines=();
    
    open(INPUT, $input_file_name) || die "Can't open $ARGV[0]\n";
    open(DUPS, ">$dup_output_file_name") || die "Error 1\n";
    open(UNIQUE, ">$unique_output_file_name") || die "Error 2\n";

        while(my $line=<INPUT>){
        chomp $line;
        my @liner=split(/\t|\s+/, $line);
   
        push @{ $lines{$liner[eval{$column-1}]}}, $line;
   
    
    
        }

        for my $key (keys %lines){
    
            my @instances=@{$lines{$key}};
                    
            if(scalar @instances==1){
        
                print UNIQUE $instances[0]."\n";
                }
            else{
                foreach(@instances){
            
             print DUPS $_."\n";
                }
              }
    
    
                }
    
    close INPUT;
    close DUPS;
    close UNIQUE;
}
sub getEgpFromCountsOutput{
    open(COUNTS_OUTPUT, "<$_[0]") || die "Can't open counts_output for parsing\n";
    open(OUTPUT, ">$_[1]");
    my $col_to_append=$_[2];
    
    while(my $line=<COUNTS_OUTPUT>){
        chomp $line;
        my @liner=split(/\t+|\s+/, $line);
      my $egp_id=shift @liner;
      my $chr=shift @liner;
        
        foreach(@liner){
            
            if ($_ ne 'NA'){
                
                print OUTPUT "$egp_id"."*".$chr."_".$_."\t$col_to_append"."\n";
                
            }
            else{
                
                last;
            }
            
            
        }
        
    }
    
}

###STEP2 SUBROUTINES
sub runStep2{
    print LOG "Running step 2\n";
        @timeData = localtime(time);
print LOG "Step 2 started: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";

    $wig_file=$_[0];
    $peak_file=$_[1];
    $sample_name=$_[2];
    
   system("cat $db_enh $peak_file |sortBed -i stdin | mergeBed -i stdin > all_merged_enhancers.txt")==0 or errorHandler($!,"Catting peaks 2 failed");
   `mkdir "wig"`;
   `mkdir "wig/$sample_name"`;
split_all_chr_wig_into_single_chr($wig_file, $sample_name, "wig/$sample_name");

#2. Create params file for Heatmap_MinMedianMax (only for the samples in the db)
`cp -r $scripts/heatmap_script/ $working_dir`;
create_params_file("$hm_dir"."/all_db_samples_params.txt", $chr_for_hm_script,"../../all_merged_enhancers.txt", "$db_wig_dir", $samples_comma_delim_string);

#3. Run Heatmap_MinMedianMax (WIG version) on all samples from db , using the master peaks file from previous step
run_hm_minmedianmax_script_new("all_db_samples_params.txt");

#4 Create params file for Heatmap_MinMedianMax  (WIG version) for the new user sample
create_params_file("$hm_dir"."/new_sample_params.txt", $chr_for_hm_script,"$working_dir/all_merged_enhancers.txt", "$working_dir"."/wig/", "$sample_name");

#5 Run Heatmap_MinMedianMax on the 10th sample, using the master peaks file
run_hm_minmedianmax_script_new("new_sample_params.txt");

#6 Column-bind the Heatmap_MinMedianMax_results for all samples
my $paste_string=" paste ";
my $add_header_command="echo -e \"chr\\tstart\\tstop\\t";
my $add_header_command2="echo -e \"enh_chr_coord\\tchr\\tstart\\tstop\\t";
my $se_header="\\t";
foreach(@samples_arr){
   $paste_string=$paste_string."$working_dir/$hm_dir/$hm_subdir/$_"."_OUTPUT/OUTPUT_ALL_CHR.txt ";
   $add_header_command=$add_header_command."$_"."_H3K4me1\\t";
   $add_header_command2=$add_header_command2."$_"."_H3K4me1\\t";
   $se_header=$se_header."$_"."_SE\\t";

}
$paste_string=$paste_string."$working_dir/$hm_dir/$hm_subdir/$sample_name"."_OUTPUT/OUTPUT_ALL_CHR.txt > all_samples_hm_min_median_max_output.txt";
system($paste_string)==0 or errorHandler($!, "Unable to cbind Heatmap output");


#5. Go through each numeric value, replace NA's with 0's and add 1 to each value
$add_header_command=$add_header_command."$sample_name"."_H3K4me1\" > ENH_floored.txt";
$add_header_command2=$add_header_command2."$sample_name"."_H3K4me1$se_header$sample_name"."_SE\"";

system("$add_header_command")==0 or errorHandler($!, "Can't append header to all_samples_hm_min_median_max_output.txt");
my $cat_str= "cat all_samples_hm_min_median_max_output.txt | sed 's/NA/0/g' | awk '{OFS=\"\t\"}{print\$1,\$2,\$3,";

for (my $i=6;$i<6*(scalar @samples_arr +1);$i+=6){
   
   $cat_str=$cat_str."\$$i+1,";
}
$cat_str=$cat_str."\$".eval{(scalar @samples_arr +1)*6}."+1}' >>ENH_floored.txt" ;
system($cat_str)==0 or errorHandler($!, "Error flooring heatmap output\n");

#6. Q-norm
print LOG "Q-norming...";
system("R --vanilla --args \"ENH_floored.txt\" \"ENH_qnormed.txt\" 12 < $scripts/qNorm.R >R_code.txt")==0 or errorHandler($!, "Unable to qNorm floored enhancers file");
print LOG "Done\n";
#7 Shannon Entropy scoring
shannon_entropy_scoring("ENH_qnormed.txt","ENH_ALL_SE_scored.txt");
 system("awk ' {printf  \"%s_%d\\n\", \$1, (\$2+\$3)/2}' ENH_ALL_SE_scored.txt | paste - ENH_ALL_SE_scored.txt  | awk '\$2 != \"chr\"' >> ENH_ALL_SE_scored2.txt " )==0 or errorHandler($!,"Error appending row names to ENH_ALL_SE_scored.txt");
system("$add_header_command2 > ENH_ALL_SE_scored.txt");
system(" cat ENH_ALL_SE_scored2.txt >> ENH_ALL_SE_scored.txt " )==0 or  errorHanlder($!,"Error renaming file");
unlink "ENH_ALL_SE_scored2.txt";
    
      @timeData = localtime(time);
print LOG "Step 2 finished: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";
  
}
###HELPER SUBROUTINES####
sub getWigFileInfo{
     my $wig_file=$_[0];
     my $wig_type="";
     my $step="";
     my $start="";
     my $chrom="";
     
     my $wig_header_lines=`grep -E 'fixedStep|variableStep' $wig_file`;
     my @wig_header_line_arr=split("\n", $wig_header_lines);
    
     foreach(@wig_header_line_arr){
          my $wig_header_line=$_;
           if($wig_header_line=~m/chrom=(\S+)\s*/){ #get chromosome info
               $chrom=$1;
           }
     
     if($wig_header_line=~m/fixedStep/){
          
         $wig_type="fixedStep";
         if($wig_header_line=~m/step=(\d+)\s*/){
          
          $step=$1;
          
         }
         if($wig_header_line=~m/start=(\d+)\s*/){
          
          $start=$1;
         }
         
         
         
     }
     elsif($wig_header_line=~m/variableStep/){
          
          $wig_type="variableStep";
          $step="NA";
          $start="NA";
     }
     else{
          
      
     undef $_;
     }
     
     if($wig_type eq "" and $step eq "" and $start eq "" and $chrom eq ""){
          undef $_; 
     }
     else{
          if (defined($_)){
               push @valid_wigs, $chrom;
          #   print "Valid chrom wig: $chrom\n"; 
               
          } 
          
     }
     
          
     }
     
     return @valid_wigs;
}
sub split_all_chr_wig_into_single_chr{
    
    my $wig_file=$_[0];
    my $sample_name=$_[1];
    my $output_dir=$_[2];
    
    open(WIG_FILE, "<$wig_file") || die "Can't open wig file: $wig_file\n";
my @chroms=getWigFileInfo($wig_file);

my %chroms_hash=();
foreach(my $i=0; $i<=$#chroms;$i++){
    local *FILE;
    open(FILE, ">$output_dir/$sample_name"."_$chroms[$i].wig") || die ("Can't open $output_dir/$sample_name"."_$chroms[$i].wig") ;
     $chroms_hash{$chroms[$i]}=*FILE;
}
my $temp_file_handle;
my $line_count=0;
while(my $line=<WIG_FILE>){
     
     if ($line=~m/chrom=(\S+)\s/){
          
          if(defined($chroms_hash{$1})){
               
               $temp_file_handle=$chroms_hash{$1};
               
               print $temp_file_handle $line;
          }
          
     }
     elsif($line =~m/^\d+/ and defined($temp_file_handle)){ # if data line
         
          print $temp_file_handle  $line;
          
     }
     
     
     
}
    
}

sub run_hm_minmedianmax_script_new{
    print LOG "Running HeatmapMinMedianMax on $_[0]\n";
   my  $params_file=$_[0];
   my $sample_name=$_[1];
    
    chdir $hm_dir;
    chdir $hm_subdir;

   system("perl getMinMedianMax.pl ../$params_file $num_cores >> log.txt ")==0 or errorHandler($!, "Heatmap script failed ");
    
    chdir $working_dir;
print LOG "\nDone\n";
}

sub create_params_file{
    print LOG "Creating params file $_[0]..."  ;
    $params_f_name=$_[0];
    my $chr_str=$_[1];
    my $abs_path_to_p_f=$_[2];
    my $abs_path_to_s_ds=$_[3];
    my $sample_dirs_str=$_[4];
   

    open(PARAMS_FILE, ">$params_f_name") || die "Can't create $params_f_name $!\n";
    
    print PARAMS_FILE "chromosomes:$chr_str\n";
    print PARAMS_FILE "absolute_path_to_peak_file:$abs_path_to_p_f\n";
    print PARAMS_FILE "absolute_path_to_sample_dirs:$abs_path_to_s_ds\n";
    print PARAMS_FILE "sample_dirs:$sample_dirs_str\n";
    
    print LOG "Done\n";
}



#STEPS 3 AND 4 SUBROUTINES##
sub runSteps3and4{
    print LOG "Running steps 3 and 4\n";
        @timeData = localtime(time);
print LOG "Steps 3 and 4 started: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";


	$user_fpkm_file=$_[0];
	$db_fpkm_file=$_[1];
	$sample_name=$_[2];
        
       
#1. Merge user genes' FPKM expression data to FPKM_all_db_lines
print LOG "Merging user genes' FPKM expression data to FPKM_all_db_lines...";
my $head1= `head -1 $db_fpkm_file`;
chomp $head1;
my $temp_header="$head1\t$sample_name"."_expr\n";
my @temp_header_arr=split(/\t+|\s+/, $temp_header);
my $new_temp_header_str=$temp_header_arr[0];
for(my $i=1;$i<$#temp_header_arr;$i++){$new_temp_header_str=$new_temp_header_str."\\t".$temp_header_arr[$i];}
$new_temp_header_str=$new_temp_header_str."\\t".$temp_header_arr[$#temp_header_arr];
system("echo -e \"$new_temp_header_str\" > FPKM_all_lines.txt");

my $merge_string="";
####if -gc not selected
if($options{'gl'}){
   system ("cat $user_fpkm_file | awk ' BEGIN { IFS=\"[ \t/\/\s]+\";OFS=\"\t\" } {print \$1\";\"\$3, \$2}' > user_fpkm_with_gene_cluster_id.txt" )==0 or errorHandler($!, "Unable to reformat user FPKM file for merging to db FPKM data");
$merge_string=" $scripts/./merge.plx -a -i $db_fpkm_file:3 user_fpkm_with_gene_cluster_id.txt:1  | awk 'BEGIN { IFS=\"[/\t/\/\s]+\";OFS=\"\t\" }  {print \$1, \$2,\$3,";
}
else{
   $merge_string=" $scripts/./merge.plx -a -i $db_fpkm_file:1 $user_fpkm_file:1  | awk ' BEGIN { IFS=\"[/\t/\/\s]+\";OFS=\"\t\" }  {print \$1, \$2,\$3,";
   
}

for(my $i=4;$i<=scalar @samples_arr +3;$i++){
   
   $merge_string=$merge_string."\$$i,"
}
$merge_string=$merge_string."\$".eval{scalar @samples_arr+5}."\}'>>FPKM_all_lines.txt";

system($merge_string)== 0 or errorHandler($!, "Unable to merge FPKM data $!\n");
print LOG "Done\n";
#2. Q-norm 
print LOG "Qnorming...";
system("R --vanilla --args \"FPKM_all_lines.txt\" \"FPKM_all_lines_qnormed.txt\" 12 < $scripts/qNorm.R >> R_code.txt ")==0 or errorHandler($!,"Unable to qNorm FPKM file");
print LOG "Done\n";
#3 Shannon Entropy scoring
shannon_entropy_scoring("FPKM_all_lines_qnormed.txt","EXPR_ALL_SE_scored.txt");
#4 merge cell-type specific genes to enhancers
print LOG "Merging cell-type specific genes to enhancer-gene pairs...";
system("echo -e \"gene_id\\tgene_coord\\tgene_name\\tgene_locus\\tenh_chr_coord\\tctcf\\tdist\"> EGP_MERGEme.txt")==0 or errorHandler($!,"Error appending header to EGP file");
system("cat d3_CTCF_all_lab.txt | sed 's/[;*]/\t/g'| awk '{printf \"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\", \$2\";\"\$3, \$1, \$2, \$3, \$4, \$5, \$6}'  | sort -u >> EGP_MERGEme.txt")==0 or errorHandler($!,"Error replacing special characters with spaces");
system("$scripts/./merge.plx -a EGP_MERGEme.txt:1 EXPR_ALL_SE_scored.txt:3 > EXP_to_pairs.txt")==0 or errorHandler($!, "Error merging expression data to EGP");
system("$scripts/./merge.plx -a EXP_to_pairs.txt:5 ENH_ALL_SE_scored.txt:1 > MERGED_prestige_FINAL.txt")==0 or errorHandler($!, "Error merging intermediate pairs data to ehnacer expression");
system("grep 'gene_id' MERGED_prestige_FINAL.txt > MERGED_prestige_FINAL2.txt");
system("grep -v 'gene_id' MERGED_prestige_FINAL.txt >> MERGED_prestige_FINAL2.txt");
system("mv MERGED_prestige_FINAL2.txt MERGED_prestige_FINAL.txt");
print LOG "Done\n";
#5 make predicted pairs files for each of cell line
print LOG "Making predicted pairs files for each cell line...";
getPredictedPairs( "MERGED_prestige_FINAL.txt",$sample_name);
print LOG "Done\n";
    @timeData = localtime(time);
print LOG "Steps 3 and 4 finished: ".eval{$timeData[4]+1}."/$timeData[3]/".eval{1900+$timeData[5]}."-$timeData[2]:$timeData[1]:$timeData[0]\n";
exit();	
}
####HELPER SUBROUTINES#######
sub getPredictedPairs{
	my $input_file_name=$_[0] ;
	my $sample_name=$_[1];
	open(INPUT, "<$input_file_name") || errorHandler($!,"Can't open master predicted pairs file");
	
my %out_files_map=();
foreach(@samples_arr){ ##creates empty output predicted pairs files for all cell lines, including the new cell line
   local *FILE;
  open(FILE, ">$_"."_ls_predicted_pairs.txt") || errorHandler ($!,"Can't create $_"."_ls_predicted_pairs.txt");
   $out_files_map{$_}=*FILE;
   #$temp_file_handle=$chroms{$chrom};
   
}
local *FILE;
open(FILE, ">$sample_name"."_ls_predicted_pairs.txt") || errorHandler ($!,"Can't create $sample_name"."_ls_predicted_pairs.txt");
$out_files_map{$sample_name}=*FILE;

my %labels_to_indeces=();
my $line_count=0;
while(my $line=<INPUT>){
   
   $line_count++;
chomp $line;
my @liner=split(/\t|\s+/, $line);

if($line_count == 1){ #if first line, get headers
   for(my $i=0;$i<scalar @liner;$i++){
      if($liner[$i] =~m/_SE/){
         push @{$labels_to_indeces{$liner[$i]}}, $i;
      }
      else{
       $labels_to_indeces{$liner[$i]}=$i;
      }
   }
   
   for my $key  ( keys %out_files_map){
      my $temp_file_handle=$out_files_map{$key};
      print $temp_file_handle  "sample_name\tgene_name\tenh_chr\tenh_start\tenh_stop\tH3K4me1\tH3K4me1_spec\tgene_expr\tgene_spec\twithin_ctcf\twithin_100KB\n";
   }
#print Dumper(%labels_to_indeces);
   
}
else{
  
  my $ctcf=$liner[$labels_to_indeces{'ctcf'}];
  my $dist=$liner[$labels_to_indeces{'dist'}];
$ctcf=~ s/g40/ctcf/;
$dist=~s/d3/<100kb/;


foreach(@samples_arr){
   my $s_name=$_;
my   $file_handle=$out_files_map{$s_name};

 print $file_handle "$s_name\t$liner[$labels_to_indeces{'gene_name'}]\t$liner[$labels_to_indeces{'chr'}]\t$liner[$labels_to_indeces{'start'}]\t$liner[$labels_to_indeces{'stop'}]\t$liner[$labels_to_indeces{$s_name.'_H3K4me1'}]\t$liner[$labels_to_indeces{$s_name.'_SE'}[1]]\t$liner[$labels_to_indeces{$s_name.'_expr'}]\t$liner[$labels_to_indeces{$s_name.'_SE'}[0]]\t$ctcf\t$dist\n" if $liner[$labels_to_indeces{$s_name.'_SE'}[0]] < 6.8 and $liner[$labels_to_indeces{$s_name.'_SE'}[1]] < 6.1 and $liner[$labels_to_indeces{$s_name.'_H3K4me1'}] > 10;

   
}

 $s_name=$sample_name;
$file_handle=$out_files_map{$s_name};
print $file_handle "$s_name\t$liner[$labels_to_indeces{'gene_name'}]\t$liner[$labels_to_indeces{'chr'}]\t$liner[$labels_to_indeces{'start'}]\t$liner[$labels_to_indeces{'stop'}]\t$liner[$labels_to_indeces{$s_name.'_H3K4me1'}]\t$liner[$labels_to_indeces{$s_name.'_SE'}[1]]\t$liner[$labels_to_indeces{$s_name.'_expr'}]\t$liner[$labels_to_indeces{$s_name.'_SE'}[0]]\t$ctcf\t$dist\n" if $liner[$labels_to_indeces{$s_name.'_SE'}[0]] < 6.8 and $liner[$labels_to_indeces{$s_name.'_SE'}[1]] < 6.1 and $liner[$labels_to_indeces{$s_name.'_H3K4me1'}] > 10;

   

}



 
   
}


my $cat_str="cat ";
foreach(@samples_arr){
   
   $cat_str="$cat_str$_"."_ls_predicted_pairs.txt ";
}
$cat_str=$cat_str." > all_db_PreSTIGE_cell_lines_ls_predicted_pairs.txt";
system("$cat_str")==0  || errorHandler($!, "Catting predicted pairs files failed");

if(defined($options{hs})){
   print LOG "Getting high stringency peaks...";
   system("awk '\$7<6 && \$9<6' all_db_PreSTIGE_cell_lines_ls_predicted_pairs.txt > all_db_PreSTIGE_cell_lines_hs_predicted_pairs.txt")==0 or errorHandler($!,"Error getting high stringency predictions for PreSTIGE cell lines");
   system("awk '\$7<6 && \$9<6' $sample_name"."_ls_predicted_pairs.txt > $sample_name"."_hs_predicted_pairs.txt")==0 or errorHandler($!,"Error getting high stringency predictions for new cell line");

}

if(defined($options{acpeaks})){
###overlap all_db_samples_predicted_pairs with acetyl peaks
print  LOG "overlapping $ac_peaks_file with predicted pairs files...";
open(ALL_DB_PREDICTED_PAIRS, "<all_db_PreSTIGE_cell_lines_ls_predicted_pairs.txt");
open(ALL_DB_PREDICTED_PAIRS_FOR_OVERLAP, ">all_db_PreSTIGE_cell_lines_ls_predicted_pairs_for_overlaps.txt");
while(my $line=<ALL_DB_PREDICTED_PAIRS>){
   chomp $line;
   if ($line !~ "within_ctcf"){
   my @liner=split(/\t|\s+/, $line);
   print ALL_DB_PREDICTED_PAIRS_FOR_OVERLAP "$liner[2]\t$liner[3]\t$liner[4]\t$liner[0]\t$liner[1]\t";
  for(my $i=5; $i<scalar @liner; $i++){print ALL_DB_PREDICTED_PAIRS_FOR_OVERLAP "$liner[$i]\t";} print ALL_DB_PREDICTED_PAIRS_FOR_OVERLAP"\n";
   }
}
my $header_line=`head -1 all_db_PreSTIGE_cell_lines_ls_predicted_pairs.txt `;
my @header=split(/\s+|\t/, $header_line);
open(TEMP, ">all_db_PreSTIGE_cell_lines_ls_predicted_pairs_overlap_with_chromatin_mark_of_interest.txt ");
print TEMP "$header[2]\t$header[3]\t$header[4]\t$header[0]\t$header[1]\t"; for(my $i=5; $i<scalar @header; $i++){print TEMP "$header[$i]\t";} print TEMP "\n";

system("intersectBed -a all_db_PreSTIGE_cell_lines_ls_predicted_pairs_for_overlaps.txt -b $ac_peaks_file -u  >> all_db_PreSTIGE_cell_lines_ls_predicted_pairs_overlap_with_chromatin_mark_of_interest.txt  ")==0 or errorHandler($!,"Unable to overlap $ac_peaks_file with $sample_name"."_ls_predicted_pairs.txt");


###overlap new_sample_predicted_pairs with acetyl peaks
open(SAMPLE_PREDICTED_PAIRS, "<$sample_name"."_ls_predicted_pairs.txt");
open(SAMPLE_PREDICTED_PAIRS_FOR_OVERLAP, ">$sample_name"."_ls_predicted_pairs_for_overlaps.txt");
while(my $line=<SAMPLE_PREDICTED_PAIRS>){
   chomp $line;
   if ($line !~ "within_ctcf"){
   my @liner=split(/\t|\s+/, $line);
   print SAMPLE_PREDICTED_PAIRS_FOR_OVERLAP "$liner[2]\t$liner[3]\t$liner[4]\t$liner[0]\t$liner[1]\t";
  for(my $i=5; $i<scalar @liner; $i++){print SAMPLE_PREDICTED_PAIRS_FOR_OVERLAP "$liner[$i]\t";} print SAMPLE_PREDICTED_PAIRS_FOR_OVERLAP"\n";
   }
}
my $header_line=`head -1 $sample_name"_ls_predicted_pairs.txt" `;
my @header=split(/\s+|\t/, $header_line);
open(TEMP, ">$sample_name"."_ls_predicted_pairs_overlap_with_chromatin_mark_of_interest.txt ");
print TEMP "$header[2]\t$header[3]\t$header[4]\t$header[0]\t$header[1]\t"; for(my $i=5; $i<scalar @header; $i++){print TEMP "$header[$i]\t";} print TEMP "\n";

system("intersectBed -a $sample_name"."_ls_predicted_pairs_for_overlaps.txt -b $ac_peaks_file -u  >> $sample_name"."_ls_predicted_pairs_overlap_with_chromatin_mark_of_interest.txt  ")==0 or errorHandler($!,"Unable to overlap $ac_peaks_file with $sample_name"."_ls_predicted_pairs.txt");


}

}

sub checkUserFPKMFileFormat{
	my $file=$_[0];
	open(FILE, "<$file");
        open(OUTPUT, ">$working_dir/user_fpkm_file.txt");
	
	while(my $line=<FILE>){
		chomp $line;
               
		if($line!~m/(.+)\t([\d]+(\.[\d]+)*)(\t.+)*/){
	errorHandler("","User FPKM file line in wrong format: $line");
	
}
                else{
                    if(! exists($options{'m'})){
                        my $val=$2;
                        $val=0 if $val<0.3;
                    print OUTPUT "$1\t".eval{$val+0.3};
                    if($options{'gl'}){
                     
                     print OUTPUT "$4";
                    }
                               
                                 print OUTPUT "\n";
                                
                    
                    }
                    else{
                     
                     print OUTPUT $line."\n";
                     
                    }
                    
                    
                }

		
	}


	return "$working_dir/user_fpkm_file.txt";
	
}

sub checkUserPeakFileFormat{
   print LOG "Checking user peak file format...";

	my $file=$_[0];
	open(FILE, "<$file");
       
	
	while(my $line=<FILE>){
		chomp $line;
		if($line!~m/^chr[^\t\s]+\t\d+(\.[\d]+)*\t\d+(\.[\d]+)*$/){
	errorHandler("","New cell line peak file has a malformed line:$line. It must be in this format: chromosome<TAB or space>peak_start<TAB or space>peak_stop");
	
}
           

		
	}


	print LOG "Done\n"
	
	
}

sub checkIfAllChromsInUserFiles{
my $peak_file=$_[0];
my $wig_file=$_[1];
my $db_enh=$_[2];
my $new_db_enh="PreSTIGE_peak_list.txt";
my $new_peak_file="new_sample_peak_list.txt";
`rm $new_db_enh` if -e $new_db_enh;
my @chr=split(/,/,$chr_for_hm_script);
my @new_chr_arr=();
foreach(@chr){
my $chr=$_;
my $in_sample_wig=`grep -w chr$chr $wig_file | wc | awk '{print \$1}'`;
chomp $in_sample_wig;
my $in_sample_peaks=`grep -w chr$chr $peak_file | wc | awk '{print \$1}'`;
chomp $in_sample_peaks;
if($in_sample_wig>0 and $in_sample_peaks>0){ push @new_chr_arr, $chr; `grep -w chr$chr $db_enh >> $new_db_enh`; `grep -w chr$chr $peak_file >> $new_peak_file`;}
if($in_sample_wig ==0 and $in_sample_peaks>0){print LOG "WARNING: chr$chr is missing from user wig file $wig_file and will be skipped\n";}
if($in_sample_wig>0 and $in_sample_peaks==0){print LOG "WARNING: chr$chr is missing from user peak file $peak_file and will be skipped\n";}
}
my $new_chr_for_hm_script= join(",", @new_chr_arr);
return ($new_chr_for_hm_script, $new_db_enh, $new_peak_file);

}


   