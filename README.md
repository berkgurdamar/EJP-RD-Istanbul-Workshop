# EJP RD Istanbul Workshop

This repository contains all required informations for performing Whole Exome Sequencing (WES) analysis.

## Steps

- Quality Control with [FastQC](https://github.com/s-andrews/FastQC)
- Trimming with [Trim-Galore](https://github.com/FelixKrueger/TrimGalore)
- Mapping with [BWA-MEM](https://github.com/lh3/bwa)
  - If necessary, merge BAMs
- Processing using [GATK](https://gatk.broadinstitute.org/hc/en-us) best practices
  - MarkDuplicates
  - FixMateInformation
  - BaseRecalibrator
  - ApplyBQSR
  - HaplotypeCaller
    - If multiple samples, CombineGVCFs
  - GenotypeGVCFs
- Annotation
  - [vcfanno](https://github.com/brentp/vcfanno)
  - [vep](https://github.com/Ensembl/ensembl-vep)
  - [vcf2db](https://github.com/quinlan-lab/vcf2db)
  - [gemini](https://github.com/arq5x/gemini)
  - [Intervar](https://github.com/WGLab/InterVar)
- Filtration, Variant Prioritization and Further Analysis

## Starting

FASTQ files are located in the '/home/projects/ejprd_istanbul_workshop/FASTQ/'.
