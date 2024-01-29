# DAPseq_pipeline_nf
Workflow written in nextflow


Input: 
1. IP sample reads (.bam); Input sample reads (.bam) 

## Pipeline Summary
1. Peak calling (MACS3)
2. Motif analysis (MEME Suite)
3. Peak annotation (HOMER)

## Requirements
```
conda environment: DAPseq_env
meme: /project/zhuzhuzhang/lyang/software/meme
```
