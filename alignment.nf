// Define parameters
params.outdiralign = "/Users/naouel/Documents/Nextflow_Rna_Seq/align_file"
params.bowtie2_exec = "/usr/local/bin/bowtie2"
params.index_path = "/Users/naouel/Documents/Nextflow_Rna_Seq/Bowtie_index/"
params.samtools_path = "/opt/anaconda3/envs/myenv2/bin/samtools"
params.threads = 4
params.srr_numbers = params.srr_numbers ?:

params.outdirtrimm = "/Users/naouel/Documents/Nextflow_Rna_seq/trimmed_files/${params.srr_numbers}_paired_trimmed.fastq.gz"

// Create channel from file pairs
reads_ch = Channel.fromFilePairs(params.outdirtrimm, checkIfExists: true)

// Define a process to handle Bowtie2 alignment
process bowtieAlignment {
    publishDir path: "${params.outdiralign}", mode: 'copy'
    input:
    tuple val(srr_numbers),path(reads)

    output:
    path "${srr_numbers}.sam"

    script:

    def idx = params.index_path + "index"
    """
    ${params.bowtie2_exec} \\
        -x "${idx}" \\
        -U "${reads}"\\
        -S "${srr_numbers}.sam"
    """
}




// Convert SAM to BAM after the alignment is completed
process samtobam {
    publishDir path: "${params.outdiralign}/${params.srr_numbers}", mode: 'copy'
    input:
    file sam

    output:
    file "output.bam"
    script:
    """
    ${params.samtools_path} view -bS $sam > "output.bam"
    """
}

// Define the main workflow
workflow {
    // Run Bowtie2 alignment
    bowtie2_out_ch = bowtieAlignment(reads_ch)

    // Convert SAM to BAM after the alignment is completed
    samtobam(bowtie2_out_ch)
}




