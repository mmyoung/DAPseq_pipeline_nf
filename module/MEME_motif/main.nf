process MEME_MOTIF {
    tag "MEME motif analysis for ${sample_id}"
    label 'motif_analysis'
    publishDir "${params.output_dir}/meme_motif", mode: 'copy'

    // conda "/project/zhuzhuzhang/lyang/software/miniconda3/envs/DAPseq_env"

    input:
    tuple val(sample_id), path(peak, stageAs: 'peaks_meme.narrowPeak')
    path fasta, stageAs: 'genome_meme.fa'

    output:
    tuple val(sample_id), path("*_meme"), optional: true, emit: meme_dir
    tuple val(sample_id), path("*.peak.fasta"), optional: true, emit: peak_fasta

    script:
    """
    module load samtools
    if [ -s peaks_meme.narrowPeak ]; then
        peak_count=\$(wc -l < peaks_meme.narrowPeak)
        echo "Processing peak file: peaks_meme.narrowPeak"
        echo "Peak file size: \${peak_count} lines"
        
        # Check if we have enough peaks for motif analysis (minimum 10 recommended)
        if [ \${peak_count} -lt 10 ]; then
            echo "Warning: Only \${peak_count} peaks found. MEME requires at least 10 peaks for reliable motif discovery."
            echo "Skipping motif analysis due to insufficient peaks."
            mkdir -p ${sample_id}_meme
            touch ${sample_id}.peak.fasta
            echo "Insufficient peaks (\${peak_count} < 10) for motif analysis" > ${sample_id}_meme/motif_analysis_skipped.txt
            exit 0
        fi
        
        # Generate fasta index if it doesn't exist
        if [ ! -f genome_meme.fa.fai ]; then
            echo "Generating FASTA index..."
            samtools faidx genome_meme.fa
        fi
        
        # Use all peaks if less than 100, otherwise take top 100
        num_peaks=\$([ \${peak_count} -lt 100 ] && echo \${peak_count} || echo 100)
        echo "Extracting sequences from top \${num_peaks} peaks..."
        
        bedtools getfasta \\
            -fi genome_meme.fa \\
            -bed <(sort -k 7 -n -r peaks_meme.narrowPeak | head -\${num_peaks} | awk 'BEGIN{OFS="\\t"}{print \$1,(\$2+\$10-30 >0)?\$2+\$10-30:1,\$2+\$10+30}') \\
            -fo ${sample_id}.peak.fasta
        
        # Check if we got sequences
        if [ -s ${sample_id}.peak.fasta ]; then
            seq_count=\$(grep -c "^>" ${sample_id}.peak.fasta)
            echo "Extracted \${seq_count} sequences for motif analysis"
            
            # Remove existing output directory to avoid conflicts
            rm -rf ${sample_id}_meme
            
            echo "Running MEME motif discovery..."
            timeout 30m meme ${sample_id}.peak.fasta \\
                -dna \\
                -revcomp \\
                -mod anr \\
                -nmotifs 2 \\
                -minw 8 \\
                -maxw 32 \\
                -maxsize 10000000 \\
                -time 1800 \\
                -o ${sample_id}_meme
            
            if [ \$? -eq 124 ]; then
                echo "MEME analysis timed out after 30 minutes"
                mkdir -p ${sample_id}_meme
                echo "MEME analysis timed out" > ${sample_id}_meme/motif_analysis_timeout.txt
            else
                echo "MEME analysis completed successfully"
            fi
        else
            echo "No sequences extracted from peaks, creating empty output"
            mkdir -p ${sample_id}_meme
            touch ${sample_id}.peak.fasta
        fi
    else
        echo "Peak file peaks_meme.narrowPeak is empty or does not exist"
        mkdir -p ${sample_id}_meme  
        touch ${sample_id}.peak.fasta
    fi
    """
}