##Dependencies/Tested with:
R 3.02
perl 5.8.8
Bedtools 2.17.0

###To INSTALL PRESTIGE###
##download bigWIGs from OSF.io into db_wig_files directory, place them in proper dirs, and convert bigWIGs to whole genome WIGs, then split into individual chromosome WIGs
 cd PRESTIGE/PRETIGE

###download from https://osf.io/9gh7b/ and place it into PRESTIGE/PRETIGE
txt_files.tar.gz
bigWIGs

##place bigWigs into proper directory structure
mv bigWIGs db_wig_files
cd db_wig_files

mkdir hg19 hg18 mm9
ls *.bw | while read file
do
genome=`echo $file | sed 's/\(hg18\|hg19\|mm9\)_\(.*\)\.bw/\1/g'`
sample=`echo $file | sed 's/\(hg18\|hg19\|mm9\)_\(.*\)\.bw/\2/g'`
test -d $genome/$sample  ||  mkdir -p $genome/$sample 
cp $file $genome/$sample 
done

##convert to WIGs, then split into chromosome WIGs
 for genome in hg18 hg19 mm9
 do
 cd $genome
 ls | while read sample_dir
 do
 cd $sample_dir
 bigWigToWig $sample_dir.bw $sample_dir.wig
 splitWigIntoChr.pl  $sample_dir.wig $sample_dir
 cd ../
 done
 cd ../
 done

 
 ##untar files
 tar -xvf scripts_for_analysis.tar.gz
tar -xvf txt_files.tar.gz
 

sh install_prestige.sh
###cd ../ and run pipeline with "perl PreSTIGE_pipeline_12_samples.pl"
