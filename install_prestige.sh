#!/bin/bash

PATH_TO_PRESTIGE_FILES=`pwd`

echo "PRESTIGE WIGS, TEXT FILES, AND SCRIPTS WILL BE PLACED IN $PATH_TO_PRESTIGE_FILES/PRESTIGE"

if [[  ! -r $PATH_TO_PRESTIGE_FILES || ! -w $PATH_TO_PRESTIGE_FILES || ! -x $PATH_TO_PRESTIGE_FILES ]]; then
echo "$PATH_TO_PRESTIGE_FILES must have read/write/execute access (chmod 777 $PATH_TO_PRESTIGE_FILES to set appropriate permissions)"
exit
fi

scripts_path="$PATH_TO_PRESTIGE_FILES/PRESTIGE/scripts_for_analysis"
txt_files_path="$PATH_TO_PRESTIGE_FILES/PRESTIGE/txt_files"
wig_files_path="$PATH_TO_PRESTIGE_FILES/PRESTIGE/db_wig_files"
hm_dir="$PATH_TO_PRESTIGE_FILES/PRESTIGE/heatmap_script"
cat PRESTIGE/PreSTIGE_pipeline_12_samples.pl | sed "s;full_path_to_scripts;$scripts_path;" |sed "s;full_path_to_wigs;$wig_files_path;" | sed "s;full_path_to_txt_files;$txt_files_path;"  > $PATH_TO_PRESTIGE_FILES/PRESTIGE/PreSTIGE_pipeline_12_samples.pl2
mv $PATH_TO_PRESTIGE_FILES/PRESTIGE/PreSTIGE_pipeline_12_samples.pl2 PRESTIGE/PreSTIGE_pipeline_12_samples.pl


cd PRESTIGE
echo "Unzipping PRESTIGE files in $PATH_TO_PRESTIGE_FILES/PRESTIGE..."
tar -xvf $PATH_TO_PRESTIGE_FILES/PRESTIGE/db_wig_files.tar.gz
tar -xvf $PATH_TO_PRESTIGE_FILES/PRESTIGE/txt_files.tar.gz
tar -xvf $PATH_TO_PRESTIGE_FILES/PRESTIGE/scripts_for_analysis.tar.gz
cd ../

