# DAPseq_pipeline_nf
A workflow for DAP-seq peak calling and related analysis.
This workflow could be used in line with the TrimAndMapping pipeline where .bam files are produced.

## The workflow includes the following steps:
1. Peak calling (MACS3)
2. Motif analysis (MEME Suite)
3. Peak annotation (HOMER)

## Test the workflow
```
nextflow run ./ -c ./test/params.config
```

## Parameters
```
--sample_sheet A tab-delimited file storing the samples information, with three columns: sample, ipbam, ctrlbam.
--fasta    Genome fasta file for the analyzing species.
--gtf    Genome gtf file for the analyzing species.
--output_dir    Name for directory for saving the results. Default: ./results
--data_dir  The folder where the raw .bam files are.
--gsize The size of analyzing genome.
```

## Requirements
```
conda environment: DAPseq_env (MACS3 and HOMER and MEME suite installed)
```
## Results 
```
test_out/
├── macs3_output
│   ├── sample1.annotatePeaks.txt
│   ├── sample1_peaks.narrowPeak
│   ├── sample1_peaks.xls
│   ├── sample1_summits.bed
│   ├── sample2_peaks.narrowPeak
│   ├── sample2_peaks.xls
│   └── sample2_summits.bed
└── meme_output
    ├── Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai
    ├── sample1_meme
    │   ├── logo1.eps
    │   ├── logo1.png
    │   ├── logo_rc1.eps
    │   ├── logo_rc1.png
    │   ├── meme.html
    │   ├── meme.txt
    │   └── meme.xml
    └── sample1.peak.fasta
```


