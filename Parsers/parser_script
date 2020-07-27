#!/bin/bash
annotated_vcf_gz=$1

# For making the SnpSift Data table.
bcftools view $annotated_vcf_gz | grep -v "start_retained_variant" | $CONDA_PREFIX/share/snpsift-*/scripts/vcfEffOnePerLine.pl | SnpSift extractFields -e "NA" - "ANN[*].GENE" "ANN[*].DISTANCE" CHROM POS ID REF ALT TYPE "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" > data/snpsiftdata.csv

# For Making the Genotype Data Table.
CHROM_POS=$(printf "CHROM\\tPOS\\t");
SAMPLE_NAMES=$(bcftools query -l $annotated_vcf_gz | paste -s -d "\t" -)
echo "$CHROM_POS$SAMPLE_NAMES"> data/Genotype_Data.csv
bcftools view $annotated_vcf_gz | bcftools query -f "%CHROM\t%POS[\t%GT]\n">> data/Genotype_Data.csv

# For making the Chromosome Name Mapping table this curently contig names in both the columns the user need to change the name of the Chromosome column if required
printf "Contig\\tChromosome\n" > data/chromosome_name_mapping.csv
bcftools view -h $annotated_vcf_gz | grep "##cont"| awk -F "=|," '{print $3 "\t" $3}' >> data/chromosome_name_mapping.csv

#printf "Sample-ID\\tPlant-ID\\tBranch-ID\\tVariety\\tGeneration\\tTreatment\\tDose\n" > data/passport.csv
#bcftools query -l $annotated_vcf_gz >> data/passport.csv

printf "Sample-ID\\tPlant-ID\\tBranch-ID\\tVariety\\tGeneration\\tTreatment\\tDose\n" > data/passport.csv
a=$(bcftools query -l $annotated_vcf_gz)
b="\tNA\tNA\tNA\tNA\tNA\tNA"
for i in ${a[*]}; do
    echo -e $i$b >> data/passport.csv;
done
