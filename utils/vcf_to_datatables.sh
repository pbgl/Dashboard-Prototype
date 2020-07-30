#!/bin/bash

#get the input file
annotated_vcf_gz=$1
#######################################

# Generate the Varaint annotation table (with SnpSift).
# we are excluding lines that contain the term "start_retained_variant", as SnpSift does not seem to understand this effect type.

bcftools view $annotated_vcf_gz | grep -v "start_retained_variant" | $CONDA_PREFIX/share/snpsift-*/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -e "NA" - "ANN[*].GENE" "ANN[*].DISTANCE" CHROM POS ID REF ALT TYPE "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" > data/snpsiftdata.tab

#######################################
# Generate a table with the genotype data

CHROM_POS=$(printf "CHROM\\tPOS\\t");
SAMPLE_NAMES=$(bcftools query -l $annotated_vcf_gz | paste -s -d "\t" -)
echo "$CHROM_POS$SAMPLE_NAMES"> data/genotype_data.tab
bcftools view $annotated_vcf_gz | bcftools query -f "%CHROM\t%POS[\t%GT]\n">> data/genotype_data.tab

#######################################
# Generate the Chromosome_Name to Contig_name mapping table. 

printf "Contig\\tChromosome\n" > data/chromosome_name_mapping.tab
bcftools view -h $annotated_vcf_gz | grep "##cont"| awk -F "=|," '{print $3 "\t" $3}' >> data/chromosome_name_mapping.tab

#######################################
#Generate the mock passport data table.

printf "Sample-ID\\tPlant-ID\\tBranch-ID\\tVariety\\tGeneration\\tTreatment\\tDose\n" > data/passport.tab
a=$(bcftools query -l $annotated_vcf_gz)
b="\tNA\tNA\tNA\tNA\tNA\tNA"
for i in ${a[*]}; do
    echo -e $i$b >> data/passport.tab;
done
