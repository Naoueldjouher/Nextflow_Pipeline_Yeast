#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = "/Users/naouel/Documents/Nextflow_Rna_seq"
params.trim_path = "/Users/naouel/Documents/Genome/Software/Trimmomatic"
params.trim_jar = "/Users/naouel/Documents/Genome/Software/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar"
params.java_path = "/usr/bin/java"
params.fastqc = "/usr/local/bin/fastqc"
params.srr_numbers = params.srr_numbers ?: params.srr_numbers  // Default value if not provided
// Update the raw pattern to use the dynamic SRR number
params.raw = "/Users/naouel/Documents/Nextflow_Rna_seq/downloaded_files/${params.srr_numbers}_*{1,2}.fastq.gz"
params.outdirtrimm = "/Users/naouel/Documents/Nextflow_Rna_seq/trimmed_files"

log.info """
R N A S E Q - N F   P I P E L I N E
===================================

trim_path: ${params.trim_path}
trim_jar: ${params.trim_jar}
java_path: ${params.java_path}
fastq: ${params.fastqc}

outdir : ${params.outdir}
srr_numbers: ${params.srr_numbers}
"""
.stripIndent()

// Create channel from file pairs
reads_ch = Channel.fromFilePairs(params.raw, checkIfExists: true)

process runTrimmomatic {

    // Specify the output directory and publish specific files
    publishDir "${params.outdirtrimm}", mode: 'copy', pattern: "${srr_numbers}_*.{fastq}"

    // Define input parameters
    input:
    tuple val(srr_numbers), path(reads)

    // Define output parameters
    output:
    file "${srr_numbers}_R1_paired_trimmed.fastq"
    file "${srr_numbers}_R2_paired_trimmed.fastq"
    file "${srr_numbers}_R1_unpaired_trimmed.fastq"
    file "${srr_numbers}_R2_unpaired_trimmed.fastq"


    script:
    """
    ${params.java_path} -jar ${params.trim_jar}  PE -threads 4  \
    ${reads[0]} ${reads[1]} \
    ${params.outdirtrimm}/${srr_numbers}_R1_paired_trimmed.fastq \
    ${params.outdirtrimm}/${srr_numbers}_R2_paired_trimmed.fastq \
    ${params.outdirtrimm}/${srr_numbers}_R1_unpaired_trimmed.fastq \
    ${params.outdirtrimm}/${srr_numbers}_R2_unpaired_trimmed.fastq \
    ILLUMINACLIP:${params.trim_path}/adapters/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    cd ${params.outdirtrimm}
    # Print some debugging information
    echo "Number of reads before trimming:"
    zcat ${reads[0]} | wc -l
    zcat ${reads[1]} | wc -l

    echo "Number of reads after trimming:"
    zcat ${params.outdirtrimm}/${srr_numbers}_R1_paired_trimmed.fastq.gz | wc -l
    zcat ${params.outdirtrimm}/${srr_numbers}_R2_paired_trimmed.fastq.gz | wc -l
    # Print file sizes for debugging
    ls -lh ${params.outdirtrimm}/${srr_numbers}_R*.fastq

    for file in ${srr_numbers}_*.fastq; do
         bgzip \$file
    done
     rm ${srr_numbers}_*.fastq
    """

}




workflow {
    // Print debug information for the input channel
    log.info "Raw files: ${params.raw}"
    reads_ch.view()

    // Directly invoke the trimmomatic process within the workflow
    trimmedFiles = runTrimmomatic(reads_ch)

}
































=======
>>>>>>> origin/main
