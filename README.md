User Guide of TMBstable
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

Introduction
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TMBstable is an innovative caller designed for the stable detection of variants.

Install
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The following packages and software need to be installed before running TMBstable. 

os

csv

math

pysam

numpy

random

Metrics

argparse

matplotlib

sklearn

scipy

skmultilearn

quickgt

samtools

bedtools

Usage
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
usage:    workFlow.py [-h] -b BAM_FILE -d INI_INFO_DIR -g REFGENOME -m SNP_METAMODEL -n SV_METAMODEL -w SNP_WINDOW -r INIREPEATFILEDIR -o OUTFILEDIR

optional arguments:

  -h, --help            show this help message and exit
  
  -b BAM_FILE, --bam_file BAM_FILE
                        input the bam file, eg: sample.bam
                        
  -d INI_INFO_DIR, --ini_info_dir INI_INFO_DIR
                        input the ini_info_dir, eg: ini_info/
                        
  -g REFGENOME, --refGenome REFGENOME
                        input the refGenome file, eg: refGenome/hg19.fa
                        
  -m SNP_METAMODEL, --snp_metaModel SNP_METAMODEL
                        the snp meta model file dir, eg: Nsnp_134.m
                        
  -n SV_METAMODEL, --sv_metaModel SV_METAMODEL
                        the sv meta model file dir, eg: Nsv_2567.m
                        
  -w SNP_WINDOW, --snp_window SNP_WINDOW
                        the snp window size, eg: 1000000
                        
  -r INIREPEATFILEDIR, --iniRepeatFileDir INIREPEATFILEDIR
                        the repeat file dir, eg: hg19RepeatFile/rmsk.txt
                        
  -o OUTFILEDIR, --outFileDir OUTFILEDIR
                        the out file dir, eg: myresult/TMBstable.vcf


Example
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Here is an example of running TMBstable on a Linux system:

python workFlow.py -b sample.bam -d /ini_info/ -g /refGenome/hg19.fa -m Nsnp_134.m -n Nsv_2567.m -w 1000000 -r /hg19RepeatFile/rmsk.txt -o /myresult/TMBstable.vcf

Output
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
The output of TMBstable is a VCF (Variant Call Format) file, specifically TMBstable.vcf, conforming to the VCFv4.2 standard. The VCF file format contains the following columns:

CHROM: The chromosome number where the variant is located.

POS: The position of the variant on the chromosome.

ID: A unique identifier for the variant, if available.

REF: The reference base(s) at the variant site.

ALT: The alternate base(s) observed at the variant site.

QUAL: Quality score of the variant call.

FILTER: Filter status of the variant, indicating if it passes quality thresholds.

INFO: Additional information about the variant, such as allele frequency, depth of coverage, and other annotations.

FORMAT: Format of the data in the genotype fields.

unknown: Sample-specific genotype information, detailing the genotype of the sample and additional metrics like genotype quality, depth, and allele count.

