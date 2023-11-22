#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.srr_numbers = params.srr_numbers ?: "SRR13978643"

params.tbi = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.vcf.gz.tbi"
params.dict = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.dict"
params.vcf = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.vcf.gz"
params.ref = "/Users/naouel/Documents/Nextflow_Rna_Seq/ref_info/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa"
params.fai = "${params.ref}.fai"
params.gatk_exec="java -Xmx4g -jar /Users/naouel/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"
params.javaCmd = "java -Xmx4g -jar"
params.bam="/Users/naouel/Documents/Nextflow_Rna_Seq/picard_gatk/${params.srr_numbers}_output.sorted.RD.RG.bam"
params.gatk_path="/Users/naouel/Documents/Nextflow_Rna_Seq/picard_gatk/Gatk"

// Create channel from file pairs
reads_ch = Channel.fromPath("${params.bam}", checkIfExists: true)
// BaseRecalibrator process
process baseRecalibrator {

    // Specify the output directory
    publishDir path: "${params.gatk_path}", mode: 'copy'

    input:
    path reads

    output:
    path "${params.srr_numbers}_recalibration_report1.grp"

    script:
    """
    ${params.gatk_exec} BaseRecalibrator \
        -I ${reads} \
        -R ${params.ref} \
        --known-sites ${params.vcf} \
        -O ${params.srr_numbers}_recalibration_report1.grp

    """
}



 process applyBQSR {

     // Specify the output directory
     publishDir path: "${params.gatk_path}", mode: 'copy'

     input:
     path reads
     path recalibration_report

     output:
     path "${params.srr_numbers}_recal.bam"

     script:
     """
     ${params.gatk_exec} ApplyBQSR \
         -R ${params.ref} \
         --bqsr ${recalibration_report} \
         -I ${reads} \
         -O ${params.srr_numbers}_recal.bam
         --create-output-bam-md5

     """
 }

workflow {
    // Run BaseRecalibrator and use its output for ApplyBQSR
    recalibrationReport = baseRecalibrator(reads_ch)
    recalibratedBAM = applyBQSR(recalibrationReport, reads_ch)
}




