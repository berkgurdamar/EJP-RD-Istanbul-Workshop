# EJP RD Istanbul Workshop

This repository contains all required information for performing Whole Exome Sequencing (WES) analysis.

## Steps

- Quality Control with [FastQC](https://github.com/s-andrews/FastQC)
- Trimming with [Trim-Galore](https://github.com/FelixKrueger/TrimGalore)
- Mapping with [BWA-MEM](https://github.com/lh3/bwa)
  - If necessary, merge BAMs
- Processing using [GATK](https://gatk.broadinstitute.org/hc/en-us) best practices
  - [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)
  - [FixMateInformation](https://gatk.broadinstitute.org/hc/en-us/articles/360036713471-FixMateInformation-Picard-)
  - [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
  - [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR)
  - [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)
    - If multiple samples, [CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037053272-CombineGVCFs)
  - [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
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

Create necessary folders

```
mkdir /home/projects/ejprd_istanbul_workshop/{user_id}/qc
mkdir /home/projects/ejprd_istanbul_workshop/{user_id}/trimming
mkdir /home/projects/ejprd_istanbul_workshop/{user_id}/mapping
mkdir /home/projects/ejprd_istanbul_workshop/{user_id}/processing
mkdir /home/projects/ejprd_istanbul_workshop/{user_id}/annotation
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
O=/home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.bam \
M=/home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup_metrics.txt \
CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/home/tmp/guest
```

#### FixMateInformation

```
picard FixMateInformation I=/home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.bam\
O=/home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.matefixed.bam \
SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/home/tmp/guest
```

#### BaseRecalibrator

```
gatk BaseRecalibrator --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.matefixed.bam \
--known-sites /home/resources/sites/hg38/dbsnp_146.hg38.vcf.gz \
--known-sites /home/resources/sites/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites /home/resources/sites/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.recal_data.table --tmp-dir /home/tmp/guest
```

#### ApplyBQSR

```
gatk ApplyBQSR --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.matefixed.bam \
--bqsr-recal-file /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.recal_data.table \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.matefixed.bqsr.bam --tmp-dir /home/tmp/guest
```

#### HaplotypeCaller

```
gatk HaplotypeCaller --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--input /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample.mdup.matefixed.bqsr.bam \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample.g.vcf.gz \
--emit-ref-confidence GVCF --annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest \
--native-pair-hmm-threads 4
```

If multiple samples

```
gatk CombineGVCFs --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample1.g.vcf.gz \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample2.g.vcf.gz \
 -O /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/combined.g.vcf.gz \
--annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest
```

#### GenotypeGVCFs

```
gatk GenotypeGVCFs --reference /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta \
--intervals /home/resources/kits/hg38/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed --interval-padding 100 \
--variant /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample.g.vcf.gz \
--output /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample_genotype.vcf.gz \
--annotation-group AS_StandardAnnotation --tmp-dir /home/tmp/guest
```


### Annotation

#### vcfanno

```
vt validate /home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample_genotype.vcf.gz \
-r /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta

vcfanno -p 4 -lua /home/projects/ejprd_istanbul_workshop/scripts/analysis_related/custom.lua \
/home/projects/ejprd_istanbul_workshop/scripts/analysis_related/config.toml \
/home/projects/ejprd_istanbul_workshop/{user_id}/variant_calling/sample_genotype.vcf.gz > /home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample_vcfanno.vcf
```

#### vep

```
/home/tools/vep/ensembl-vep-release-109/vep -i /home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample_vcfanno.vcf \
--offline --cache --dir /home/tools/vep/ --cache_version 109 \
--assembly GRCh38 --hgvsg --everything --total_length --fork 4 \
--force_overwrite -o /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample_vep.vcf \
--vcf --fasta /home/resources/reference/homo_sapiens/hg38/ucsc.hg38.fasta
```

#### vcf2db

```
vcf2db.py /home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample_vep.vcf \
/home/projects/ejprd_istanbul_workshop/sample.ped /home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.gemini.db
```

#### gemini

```
/home/tools/anaconda3/envs/wgs_gemini/bin/gemini query -q \
            	"select gene, chrom, start, end, ref, alt, type, \
            	gt_quals, gt_depths, gt_alt_depths, gt_alt_freqs, \
            	gt_types, \
            	is_exonic, is_coding, is_lof, is_splicing, \
            	hgvsc, hgvsp, codon_change, aa_change, biotype, impact_so, impact_severity, \
            	polyphen_pred, sift_pred, max_af, gnomade_af, clinvar_sig, \
            	clinvar_disease_name, \
            	existing_variation, \
            	pubmed \
            	from variants where impact_severity == 'HIGH' or impact_severity == 'MED' \
            	or polyphen_pred == 'possibly_damaging' or \
            	polyphen_pred == 'probably_damaging' or sift_pred == 'deleterious' or \
            	sift_pred == 'deleterious_low_confidence'" \
              --header /home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.gemini.db > \
              /home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.family.query.txt
```

#### InterVar

```
python /home/tools/InterVar/Intervar.py --buildver=hg38 \
        	--input=/home/projects/ejprd_istanbul_workshop/{user_id}/processing/sample_vep.vcf \
        	--database_intervar=/home/tools/InterVar/intervardb/ \
        	--table_annovar=/home/tools/annovar/table_annovar.pl \
        	--convert2annovar=/home/tools/annovar/convert2annovar.pl \
        	--annotate_variation=/home/tools/annovar/annotate_variation.pl \
        	--database_locat=/home/tools/annovar/humandb/ \
        	--input_type=VCF \
        	--output=/home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample
```

### Create Final Report

```
R CMD BATCH /home/projects/ejprd_istanbul_workshop/scripts/final_report_script.R \
/home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.hg38_multianno.txt.intervar \
/home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.query.txt \
/home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.hg38_multianno.txt \
/home/projects/ejprd_istanbul_workshop/{user_id}/annotation/sample.final.annotation.tsv
```
