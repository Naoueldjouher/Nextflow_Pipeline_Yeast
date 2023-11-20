#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.outdir = "/Users/naouel/Documents/Nextflow_Rna_seq"  // Update the output directory
params.sratoolkit = "/Users/naouel/Documents/Genome/Software/sratoolkit.3.0.7-mac64/bin"
params.srr_numbers = "/Users/naouel/Documents/Nextflow_Rna_Seq/SRR_Acc_List.txt" // Add the SRR numbers you want to download

log.info """\

        R N A S E Q - N F   P I P E L I N E
        ===================================

        reference genome : ${params.ref_genome}
        bowtie2: ${params.bowtie2}
        sratoolkit: ${params.sratoolkit}

        outdir : ${params.outdir}
        """
        .stripIndent()

/*
 * define the `download` process that downloads SRR numbers
 * using the SRA Toolkit and bgzip
 */

process DOWNLOAD {

    input:
    val srr_numbers
    path outdir

    output:
    path "${outdir}/downloaded_files/*.fastq.gz"

   script:
   """
   mkdir -p ${outdir}/downloaded_files
   cat ${srr_numbers} | while read srr_number; do
   ${params.sratoolkit}/fastq-dump --outdir ${outdir}/downloaded_files --split-files $srr_number
   cd ${outdir}/downloaded_files
   for file in *.fastq; do
   bgzip \$file
   done
   done
   """
}

workflow {

    download_ch = DOWNLOAD(params.srr_numbers, params.outdir)

}



















