# DAPseq_pipeline_nf
A workflow for DAP-seq peak calling and related analysis.
Trying to modulize the workflow of DAP analysis to save repetitive work.

## The workflow includes the following steps:
1. Reads trimming (trim_galore)
2. Clean reads mapping (bowtie2)
1. Peak calling (MACS3)
   ```
   macs3 \\
        callpeak \\
        --gsize $macs2_gsize \\
        --format $format \\
        --name $prefix \\
        --call-summits \\
        --keep-dup 1 \\
        --treatment $ipbam \\
        $CONTROL_CMD
   
   ```
3. Motif analysis (MEME Suite)
   
   Filtering the top 100 strongest (by peak score) peaks for motif finding (summit +-30bp).
   ```
   meme ${sample_id}.peak.fasta \\
        -dna \\
        -revcomp \\
        -mod anr \\
        -nmotifs 2 \\
        -minw 8 \\
        -maxw 32 \\
        -maxsize 10000000 \\
        -time 1800 \\
        -o ${sample_id}_meme
   ```
5. Peak annotation (HOMER)
   
   Assign peak to their closest genes.
   
7. FRIP score calculation (shell)
9. Produce final report (Rmarkdown)

## Prerequisite
```
nextflow program
conda environment: DAPseq_env (trim_galore and bedtools and bowtie2 and MACS3 and HOMER and MEME suite installed)
R packages: dplyr, DT, kableExtra, base64enc
Indexed genome

```

## Test the workflow
```
git clone git@github.com:mmyoung/DAPseq_pipeline_nf.git
conda activate DAPseq_env
nextflow run /project/gzy8899/lyang/DAPseq_pipeline_nf -params-file params.yml  ## yml file saving all parameters, refer to the file in ./test folder for format

## if the run was terminated in the middle without specific error, try:
nextflow run /project/gzy8899/lyang/DAPseq_pipeline_nf -params-file params.yml -resume
```

## Parameters

Parameters can be passed to the pipeline in the command line or be put in a .yml file and passed to the pipeline in a whole (-params-file params.yml)
```
--fq_sheet A tab-delimited file storing the samples information, with five columns: sample,fq1,fq2,single_end,control
--fasta    Genome fasta file for the analyzing species.
--gtf    Genome gtf file for the analyzing species.
--output_dir    Name for directory for saving the results. Default: ./results
--fq_dir  The folder where the raw .fastq files are.
--gsize The size of analyzing genome.
--bowtie_idx The bowtei2 index directory. ## built with bowtie2-build command
--gsize The size of analyzing genome. ## parameter for MACS2
--control_samples ## comma separated sample ID that are controls in the peak calling program, those samples won't be used for peak detection, if there are no contral samples, put null

```

**Caveats**:
* need to go through the scripts to make sure the path to softwares are executable for current user.
* pay attention to the path to cutadapt tool which is required to execute trim_galore and installation of it isn't compatibale with the DAPseq_env environment

## example of sample sheet passed to --fq_sheet
fq_sheet.csv (comma separated meta data)
```
sample,fq1,fq2,single_end,control
IP,SRR27496336_1.fastq,SRR27496336_2.fastq,0,Input
Input,SRR27496337_1.fastq,SRR27496337_2.fastq,0,
```
sample: sample name or id, need to be unique
fq1: read1 fastq name for pair-end sequencing (the only fastq file name for single-end sequencing)
fq2: read2 fastq name for pair-end sequencing (empty for single-end sequencing)
single_end: indicate whether the data is from single-end sequencing (1: single-end; 0: pair-end)
control; the corresponding control sample name/id for the DAP library (empty if there is no control)

## Results 
The intermediate output from all procedures are saved in the output directory the user defined, the structure looks like the following:

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
│   ├── Input_val_1.fq.gz
│   ├── Input_val_2.fq.gz
│   ├── IP_val_1.fq.gz
│   ├── IP_val_2.fq.gz
│   ├── SRR27496336_1.fastq_trimming_report.txt
│   ├── SRR27496336_2.fastq_trimming_report.txt
│   ├── SRR27496337_1.fastq_trimming_report.txt
│   └── SRR27496337_2.fastq_trimming_report.txt
└── report
    ├── read_peak.num.summary
    └── report.html
    
```
The report.html output is the one to check with some important statistics report about the data.

## Report explanation
```
sample: sample names
raw reads pairs#: total raw read (pairs for pair-end sequencing)
clean read pairs#: total read pairs passing trim_galore filters (-q 20, --length 20)
mapping ratio: % of clean reads that are mapped to genome
unique mapped pairs#: total reads number that are uniquely mapped to genome (filtered with samtools view -q 20)
peak#: total peak number output from MACS3
min5fold peak#: number of peaks that have peak score >5 in .narrowPeak output from MACS3
max peak score: the highest peak score in .narrowPeak file
peak reads#: totol number of reads that falls into peaks
FRIP_score: total percentage of reads that falls into peaks
motif: output logo from meme
```

