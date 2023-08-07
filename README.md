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

FASTQ files are located in the `/home/projects/ejprd_istanbul_workshop/FASTQ/`.

Activate the conda environment
```
conda activate EpiRARE_new
```

### Quality Control

```
fastqc --dir /home/tmp/guest \
--outdir /home/projects/ejprd_istanbul_workshop/{user_id}/qc \
--threads 4 --quiet --noextract \
/home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_1.fastq.gz /home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_2.fastq.gz
```

### Trimming

```
trim_galore --phred33 --quality 20 --gzip --length 35 \
--trim-n --output_dir /home/projects/ejprd_istanbul_workshop/{user_id}/trimming \
--retain_unpaired --cores 4 \
--paired /home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_1.fastq.gz /home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_2.fastq.gz
```

### Mapping

```
bwa mem -R "\@RG\\tID:{lane}\\tPL:ILLUMINA\\tLB:Twist_Comprehensive\\tSM:sample\" \
-M -t 8 /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
/home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_1.fastq.gz \
/home/projects/ejprd_istanbul_workshop/FASTQ/sample_L001_2.fastq.gz > /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001.sam

sambamba view --nthreads=4 --with-header --show-progress --sam-input --format=bam \
--output-filename=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001_unsorted.bam \
/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001.sam

rm -f /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001.sam

sambamba sort --tmpdir /home/tmp/guest --nthreads=4 \
--out=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.bam \
/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001_unsorted.bam

rm -f /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001_unsorted.bam

# picard BuildBamIndex INPUT=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.bam TMP_DIR=/home/tmp/guest VALIDATION_STRINGENCY=LENIENT
```

Combine multiple lanes

```
sambamba merge -t 4 /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.bam \
/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L001.bam \
/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_L002.bam

picard BuildBamIndex INPUT=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.bam \
TMP_DIR=/home/tmp/guest VALIDATION_STRINGENCY=LENIENT
```

### Processing using GATK best practices

#### MarkDuplicates

```
picard MarkDuplicates I=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.bam \
O=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.bam \
M=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup_metrics.txt \
CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/tmp/guest
```

#### FixMateInformation

```
picard FixMateInformation I=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.bam\
O=/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.matefixed.bam \
SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/tmp/guest
```

#### BaseRecalibrator

```
gatk BaseRecalibrator --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.matefixed.bam \
--known-sites /home/resources/sites/hg38/dbsnp_146.hg38.vcf.gz \
--known-sites /home/resources/sites/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites /home/resources/sites/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.recal_data.table --tmp-dir /home/tmp/guest
```

#### ApplyBQSR

```
gatk ApplyBQSR --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.matefixed.bam \
--bqsr-recal-file /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.recal_data.table \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.matefixed.bqsr.bam --tmp-dir /home/tmp/guest
```

#### HaplotypeCaller

```
gatk HaplotypeCaller --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.mdup.matefixed.bqsr.bam \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.g.vcf.gz \
--emit-ref-confidence GVCF --annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest \
--native-pair-hmm-threads 4
```

If multiple samples

```
gatk CombineGVCFs --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample1.g.vcf.gz \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample2.g.vcf.gz \
 -O /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/combined.g.vcf.gz \
--annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest
```

#### HaplotypeCaller

```
gatk GenotypeGVCFs --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample.g.vcf.gz \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_genotype.vcf.gz \
--annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest
```


### Annotation

#### vcfanno

```
vt validate /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_genotype.vcf.gz \
-r /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta

vcfanno -p 4 -lua /home/projects/ejprd_istanbul_workshop/scripts/analysis_related/custom.lua \
/home/projects/ejprd_istanbul_workshop/scripts/analysis_related/config.toml \
/home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_genotype.vcf.gz > /home/projects/ejprd_istanbul_workshop/{user_id}/mapping/sample_vcfanno.vcf
```
