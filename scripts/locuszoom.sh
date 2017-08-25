#!/bin/bash
#PBS -l walltime=00:01:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -M smozaffari@uchicago.edu
#PBS -m e
#PBS -N $name

# Katie Igartua
#this file creates Hutterite LD file for Locus zoom for SNP in question
# variables need at qsub are: imputationID, rsID, chr, gwas_file (.assoc.txt from gemma) and output_folder. 

# example -v imputationID=11448086,rsID=rs17242388,chr=19,gwas_file=/group/ober-resources/users/cigartua/Hutterites_GeneBased/2016/gemma/LDL/output/LDL.all.URE_CR_0.5.assoc.txt,pheno=LDL,output_file=/group/ober-resources/users/cigartua/Hutterites_GeneBased/2016/gemma/$pheno/output 
#make_generic_locuszoom_Sahar.sh

#Katie's script - modified for new imputation
#08.25.17

module load gcc/6.2.0 
module load R
module load plink
 
#imputationID=$1
#rsID=$2
#chr=$3
#gwas_file=$4
#pheno=$5
#output_folder=$6


cd /group/ober-resources/users/smozaffari/LocusZoom_newimputation
mkdir $imputationID
cd $imputationID

#calculate LD within 1MB of query SNP
plink --bfile /group/ober-resources/users/smozaffari/LocusZoom_newimputation/data/qc --ld-snp $imputationID -ld-window 1000000 --make-founders --out $imputationID.r2 --r2

#format Plink Output
sed 's/ \+ /\t/g'  $imputationID.r2.ld  > $imputationID.r2.ld.tab.temp

#grab colums required for LD locuszoom file
awk -F "\t" 'BEGIN {OFS="\t"} {print $2":"$3,$5":"$6,$8}' $imputationID.r2.ld.tab.temp >  $imputationID.r2.ld.tab_chrpos
#awk -F "\t" 'BEGIN {OFS="\t"} {print $4,$7,$8}' $imputationID.r2.ld.tab.temp >  $imputationID.r2.ld.tab

#convert ImputationIDs to rsIDs
awk -F "\t"  'BEGIN {OFS="\t"} NR==FNR {h[$46] = $16; next} {if(h[$2]) {print h[$1],h[$2],0,$3}}' /group/ober-resources/ $imputationID.r2.ld.tab > $imputationID.r2.ld.tab_rsID
#place header in ld file
#cat /group/ober-resources/users/smozaffari/LocusZoom_newimputation/scripts/locuszoom_header $imputationID.r2.ld.tab_rsID > $imputationID.r2.ld.tab_rsID_locuszoom_input.txt
cat /group/ober-resources/users/smozaffari/LocusZoom_newimputation/scripts/locuszoom_header $imputationID.r2.ld.tab_chrpos > $imputationID.r2.ld.tab_rsID_locuszoom_input.txt    

 
#get GWAS file ready
ss=$(echo $gwas_file | sed 's/gz/txt/g')

echo $ss

awk -v chr="$chr" '$1==chr' <(gzip -dc  $gwas_file ) >  $ss.$chr
gwas_file=$ss

#awk -v chr="$chr" '$1==chr'  $gwas_file  >  $ss.$chr
#gwas_file=$ss

 
sed 's/ /\t/g'  $gwas_file.$chr >  $gwas_file.$chr.tab
rm  $gwas_file.$chr
#for gemma files
#awk -F "\t"  'BEGIN {OFS="\t"} NR==FNR {h[$46] = $16; next} {if(h[$1]) {print h[$1], $11}}' /group/ober-resources/users/rlee/hutt_annotation/annotation_data/qc.imputed_cgi.annovar_plink_annotations.hg19_multianno.txt $gwas_file.$chr.tab > $gwas_file.$chr.locuszoom

#for POGWAS files:
#cut -f1,14 $gwas_file.$chr.tab > $gwas_file.$chr.locuszoom


cut -f2,20 $gwas_file.$chr.tab > $gwas_file.$chr.locuszoom  
head $gwas_file.$chr.locuszoom
cat  /group/ober-resources/users/smozaffari/LocusZoom/scripts/metal_header  $gwas_file.$chr.locuszoom >  $gwas_file.$chr.locuszoom.txt 

# run locuzoom 
echo "/group/ober-resources/users/smozaffari/LocusZoom/locuszoom/bin/locuszoom --metal $gwas_file.$chr.locuszoom.txt  --ld $imputationID.r2.ld.tab_rsID_locuszoom_input.txt --refsnp $rsID --build hg19 theme=publication title=$pheno ldColors=\"#F1F8FD,#F6C667,#FF7C38,#E03E36,#B80D57,#700961\" ldCuts=\"0,.2,.4,.6,.8,1\" condLdLow=\"#F6C667\"  --flank 1MB --prefix $name rfrows=10"
/group/ober-resources/users/smozaffari/LocusZoom/locuszoom/bin/locuszoom --metal $gwas_file.$chr.locuszoom.txt  --ld $imputationID.r2.ld.tab_rsID_locuszoom_input.txt --refsnp $rsID --build hg19 theme=publication title=$pheno ldColors="#F1F8FD,#F6C667,#FF7C38,#E03E36,#B80D57,#700961" ldCuts="0,.2,.4,.6,.8,1" condLdLow="#F6C667"  --flank 1MB --prefix $name rfrows=10

#cd /group/ober-resources/users/smozaffari/LocusZoom_newimputation/plots
#mv $imputationID $output_folder
