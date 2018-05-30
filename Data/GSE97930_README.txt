Details on adaptors and cell barcodes used are available from http://mccarrolllab.com/dropseq/:
- Online-Dropseq-Protocol-v.-3.1-Dec-20152.pdf

Batch script to generate UMI count matrix (cell barcodes by gene counts) from fastq files:

#!/bin/sh

dropseq_dir=
star_dir=
picard_dir=
ref_dir=
refIndex=
refFasta=
refGTF=
cutadapt_dir=
kneeplot_dir=
species=
cells=
sampleName=
R1=
R2=

$cutadapt_dir/cutadapt -a TTTTTTTTTT$ -m 20 -e 0.4 --discard-untrimmed -o $R1 -p $R2 $sampleName'_R1_001.fastq.gz' $sampleName'_R2_001.fastq.gz'

#extract cell barcodes and UMI from read1, convert pair-end fastq files to bam file, and generate fastq file for alignment
java -jar $picard_dir/picard.jar FastqToSam FASTQ=$R1 FASTQ2=$R2 SAMPLE_NAME=$sampleName OUTPUT=/dev/stdout | \
java -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=$sampleName'.unaligned.sorted.bam' SORT_ORDER=queryname
$dropseq_dir/TagBamWithReadSequenceExtended I=$sampleName'.unaligned.sorted.bam' O=/dev/stdout SUMMARY=$sampleName'.cell_tag_report.txt' BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=$sampleName'.molecule_tag_report.txt' BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/FilterBAM TAG_REJECT=XQ I=/dev/stdin O=/dev/stdout | \
$dropseq_dir/TrimStartingSequence INPUT=/dev/stdin OUTPUT=/dev/stdout OUTPUT_SUMMARY=$sampleName'.adapter_trimming_report.txt' SEQUENCE='AAGCAGTGGTATCAACGCAGAGTGAATGGG' MISMATCHES=0 NUM_BASES=5 | \
$dropseq_dir/PolyATrimmer INPUT=/dev/stdin OUTPUT=/dev/stdout MISMATCHES=0 NUM_BASES=6 | \
tee $sampleName'.unaligned.tagged.bam' | \
java -jar $picard_dir/picard.jar SamToFastq INPUT=/dev/stdin FASTQ=$sampleName'.unaligned.tagged.fastq'

#mapping
$star_dir/STAR --runThreadN 8 --genomeDir $ref_dir/$refIndex --readFilesIn $sampleName'.unaligned.tagged.fastq' --outFileNamePrefix $sampleName. --outSAMunmapped Within

java -jar $picard_dir/picard.jar SortSam I=$sampleName'.Aligned.out.sam' O=/dev/stdout SO=queryname | \
java -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir/$refFasta UNMAPPED_BAM=$sampleName'.unaligned.tagged.bam' ALIGNED_BAM=/dev/stdin O=/dev/stdout INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false| \
$dropseq_dir/TagReadWithGene I=/dev/stdin O=$sampleName'.aligned.gene.bam' ANNOTATIONS_FILE=$ref_dir/$refGTF TAG=GE
$dropseq_dir/DetectBeadSynthesisErrors I=$sampleName'.aligned.gene.bam' O=$sampleName'.aligned.clean.bam' OUTPUT_STATS=$sampleName'.synthesis_stats.txt' SUMMARY=$sampleName'.synthesis_stats_summary.txt' NUM_BARCODES=$cells PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
$dropseq_dir/BAMTagHistogram I=$sampleName'.aligned.clean.bam' O=$sampleName'.readsByBarcode.txt.gz' TAG=XC
$dropseq_dir/DigitalExpression I=$sampleName'.aligned.clean.bam' O=$sampleName'.counts.tsv' SUMMARY=$sampleName'.count_summary.txt' NUM_CORE_BARCODES=$cells EDIT_DISTANCE=1 
