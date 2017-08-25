#!/bin/bash

#PBS -l walltime=01:00:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -M cigartua@uchicago.edu
#PBS -m abe

# Katie Igarta
#this file creates Hutterite LD file for Locus zoom for SNP in question
# variables need at qsub are: imputationID, rsID, chr, gwas_file (.assoc.txt from gemma) and output_folder. 

# example -v imputationID=11448086,rsID=rs17242388,chr=19,gwas_file=/group/ober-resources/users/cigartua/Hutterites_GeneBased/2016/gemma/LDL/output/LDL.all.URE_CR_0.5.assoc.txt,pheno=LDL,output_file=/group/ober-resources/users/cigartua/Hutterites_GeneBased/2016/gemma/$pheno/output 
#make_generic_locuszoom_Sahar.sh
module load R
module load plink
 

cd /group/ober-resources/users/cigartua/MICROBIOME/Hutterites_Nasal/qiime/mapping/locuszoom
mkdir $imputationID
cd $imputationID
#calculate LD within 1MB of query SNP
plink --bfile /group/ober-resources/users/cigartua/Hutterite_annotation/qc --ld-snp $imputationID -ld-window 1000000 --make-founders --out $imputationID.r2 --r2

#format Plink Output
sed 's/ \+ /\t/g'  $imputationID.r2.ld  > $imputationID.r2.ld.tab.temp

#grab colums required for LD locuszoom file
awk -F "\t" 'BEGIN {OFS="\t"} {print $4,$7,$8}' $imputationID.r2.ld.tab.temp >  $imputationID.r2.ld.tab

#convert ImputationIDs to rsIDs
awk -F "\t"  'BEGIN {OFS="\t"} NR==FNR {h[$46] = $16; next} {if(h[$2]) {print h[$1],h[$2],0,$3}}' /group/ober-resources/users/rlee/hutt_annotation/annotation_data/qc.imputed_cgi.annovar_plink_annotations.hg19_multianno.txt $imputationID.r2.ld.tab > $imputationID.r2.ld.tab_rsID
#place header in ld file
cat ../locuszoom_header $imputationID.r2.ld.tab_rsID > $imputationID.r2.ld.tab_rsID_locuszoom_input.txt
  
#get GWAS file ready
awk -v chr="$chr" '$1==chr' $gwas_file >  $gwas_file.$chr
 
sed 's/ /\t/g'  $gwas_file.$chr >  $gwas_file.$chr.tab
rm  $gwas_file.$chr
awk -F "\t"  'BEGIN {OFS="\t"} NR==FNR {h[$46] = $16; next} {if(h[$2]) {print h[$2], $11}}' /group/ober-resources/users/rlee/hutt_annotation/annotation_data/qc.imputed_cgi.annovar_plink_annotations.hg19_multianno.txt $gwas_file.$chr.tab > $gwas_file.$chr.locuszoom
cat  /group/ober-resources/users/cigartua/MICROBIOME/Hutterites_Nasal/qiime/mapping/locuszoom/metal_header  $gwas_file.$chr.locuszoom >  $gwas_file.$chr.locuszoom.txt 

# run locuzoom 
/group/ober-resources/users/cigartua/bin/locuszoom/bin/locuszoom --metal $gwas_file.$chr.locuszoom.txt  --ld $imputationID.r2.ld.tab_rsID_locuszoom_input.txt --refsnp $rsID --build hg19 --flank 500kb

cd /group/ober-resources/users/cigartua/MICROBIOME/Hutterites_Nasal/qiime/mapping/locuszoom
mv $imputationID $output_folder
