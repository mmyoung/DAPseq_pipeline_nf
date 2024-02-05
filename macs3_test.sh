#!/bin/bash
#SBATCH --job-name=bismark_report
#SBATCH --output=bismark_report.out
#SBATCH --error=bismark_report.err
#SBATCH --account=pi-zhuzhuzhang
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=ALL  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lyang13@uchicago.edu  # Where to send email

cd /project/gzy8899/lyang/DAPseq_pipeline_nf
sample_id=test
peak=/project/gzy8899/lyang/DAPseq_pipeline_nf/test_out/macs3_output/sample1_peaks.narrowPeak
fasta=/project/gzy8899/references/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
conda init
source activate DAPseq_env
#/project/zhuzhuzhang/lyang/software/bedtools2/bin/bedtools getfasta -fi ${fasta} -bed ${peak} -fo ${sample_id}.peak.fasta
meme ${sample_id}.peak.fasta -dna -o ${sample_id}_meme
