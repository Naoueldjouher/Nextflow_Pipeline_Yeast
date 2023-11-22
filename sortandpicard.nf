#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define parameters
params.samtools_path = "/opt/anaconda3/envs/myenv2/bin/samtools"
params.threads = 4
params.srr_numbers = params.srr_numbers ?:
params.align = "/Users/naouel/Documents/Nextflow_Rna_Seq/align_file/"
params.picard = " /opt/anaconda3/envs/myenv2/share/picard-2.27.5-0/picard.jar"
params.javaCmd =  "/opt/anaconda3/envs/myenv2/bin/java -Xmx4g -jar"
params.ref = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa"
params.fai = "${params.ref}.fai"
params.tbi = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.vcf.gz.tbi"
params.dict = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.dict"
params.vcf = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.vcf.gz"
params.picard_gatk = "/Users/naouel/Documents/Nextflow_Rna_Seq/picard_gatk"

// Create channel from file pairs
reads_ch = Channel.fromPath("${params.align}/${params.srr_numbers}/output.bam", checkIfExists: true)
// Add sort process
process sortBam {
    // Use the same work directory as addReadGroupProcess
    storeDir "${params.picard_gatk}/work"

    // Specify the output directory
    publishDir path: "${params.picard_gatk}", mode: 'copy'

    // Specify the input BAM file
    input:
    file reads

    // Output file
    output:
    file "${params.srr_numbers}_output.sorted.bam"

    script:
    """
    ${params.samtools_path} sort -o ${params.srr_numbers}_output.sorted.bam ${reads}
    """
}

// Add read group process
process addReadGroupProcess {
    // Use a separate work directory for this process
    storeDir "${params.picard_gatk}/work"

    // Specify the output directory
    publishDir path: "${params.picard_gatk}", mode: 'copy'

    // Input channel
    input:
    file bamfile

    // Output file
    output:
    file "${params.srr_numbers}_output.sorted.RD.RG.bam"

    script:
    """
   echo "I=${bamfile}"
   echo "O=${params.srr_numbers}_output.sorted.RD.RG.bam"
   echo "RGLB=${params.srr_numbers}"

   ${params.javaCmd} ${params.picard} AddOrReplaceReadGroups \
   I=${bamfile} \
   O=${params.srr_numbers}_output.sorted.RD.RG.bam \
   RGID=1 RGLB=${params.srr_numbers} RGPL=ILLUMINA RGPU=run RGSM=${params.srr_numbers}

   """


}



// Add index process with dependency on sortBam
process indexBam {
    // Use the same work directory as addReadGroupProcess
    storeDir "${params.picard_gatk}/work"

    // Specify the output directory
    publishDir path: "${params.picard_gatk}", mode: 'copy'

    // Specify the input sorted BAM file
    input:
    file sortedBamFile

    // Output file
    output:
    file "${params.srr_numbers}_output.sorted.RD.RG.sorted.bam.bai"

    script:
    """
    ${params.samtools_path} index ${sortedBamFile} ${sortedBamFile}.bai
    """
}


// Run the workflow
workflow {
    // Run Picard MarkDuplicates workflow
    sorted_bam = sortBam(reads_ch)
    picard_index = addReadGroupProcess(sorted_bam)
    indexBam(picard_index)
}








