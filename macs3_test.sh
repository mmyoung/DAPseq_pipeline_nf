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
conda activate DAPseq_env

macs3 \
      callpeak \
      --gsize 1000 \
      --format BAMPE \
      --name sample2 \
      --treatment /project/gzy8899/lyang/DAPseq_pipeline_nf/test/bam_file/sample2_sorted.bam \
      --control /project/gzy8899/lyang/DAPseq_pipeline_nf/test/bam_file/sample1_sorted.bam