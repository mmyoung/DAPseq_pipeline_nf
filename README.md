# DAPseq_pipeline_nf
A workflow for DAP-seq peak calling and related analysis.
This workflow could be used in line with the TrimAndMapping pipeline where .bam files are produced.

## The workflow includes the following steps:
1. Peak calling (MACS3)
2. Motif analysis (MEME Suite)
3. Peak annotation (HOMER)

## Test the workflow
```
nextflow run ./ -params-file params.yml
```

## Parameters
```
--sample_sheet A tab-delimited file storing the samples information, with five columns: sample,fq1,fq2,single_end,control
--fasta    Genome fasta file for the analyzing species.
--gtf    Genome gtf file for the analyzing species.
--output_dir    Name for directory for saving the results. Default: ./results
--fq_dir  The folder where the raw .fastq files are.
--gsize The size of analyzing genome.
--bowtie_idx The bowtei2 index directory.
--prime5_trim_len How many bases to trim for the 5' of reads.
--prime3_trim_len How many bases to trim for the 3' of reads.
--gsize The size of analyzing genome.

```

## Requirements
```
conda environment: DAPseq_env (MACS3 and HOMER and MEME suite installed)
Indexed genome
```

## Input, example
1. fq_sheet.csv
```
sample,fq1,fq2,single_end,control
IP,SRR27496336_1.fastq,SRR27496336_2.fastq,0,Input
Input,SRR27496337_1.fastq,SRR27496337_2.fastq,0,
```

## Results 
```
├── alignment
│   ├── Input.bowtie2.log
│   ├── Input_map_sorted.bam
│   ├── IP.bowtie2.log
│   └── IP_map_sorted.bam
├── bw_output
│   ├── Input_sorted_bam.bw
│   └── IP_sorted_bam.bw
├── coverage_out
│   ├── Input_base.depth
│   ├── Input_depth.pdf
│   ├── IP_base.depth
│   └── IP_depth.pdf
├── deduplicateion_out
│   ├── Input_dedup_Q20_sorted.bam
│   ├── Input_dedup_Q20_sorted.bam.bai
│   ├── IP_dedup_Q20_sorted.bam
│   └── IP_dedup_Q20_sorted.bam.bai
├── macs3_output
│   ├── IP.annotatePeaks.txt
│   ├── IP_peaks.narrowPeak
│   ├── IP_peaks.xls
│   └── IP_summits.bed
├── meme_output
│   ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
│   ├── IP_meme
│   │   ├── logo1.eps
│   │   ├── logo1.png
│   │   ├── logo_rc1.eps
│   │   ├── logo_rc1.png
│   │   ├── meme.html
│   │   ├── meme.txt
│   │   └── meme.xml
│   └── IP.peak.fasta
└── trimm
    ├── Input_val_1.fq.gz
    ├── Input_val_2.fq.gz
    ├── IP_val_1.fq.gz
    ├── IP_val_2.fq.gz
    ├── SRR27496336_1.fastq_trimming_report.txt
    ├── SRR27496336_2.fastq_trimming_report.txt
    ├── SRR27496337_1.fastq_trimming_report.txt
    └── SRR27496337_2.fastq_trimming_report.txt
    
```


