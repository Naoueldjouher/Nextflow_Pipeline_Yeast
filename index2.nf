params.ref_genome = "/Users/naouel/Documents/Nextflow_Rna_Seq/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa"
params.bowtie2 = "/usr/local/bin/bowtie2"
params.outdir = "/Users/naouel/Documents/Nextflow_Rna_Seq/Bowtie_index"

log.info """
        R N A S E Q - N F   P I P E L I N E
        ===================================

        reference genome: ${params.ref_genome}
        bowtie2: ${params.bowtie2}

        outdir: ${params.outdir}
        """.stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the reference genome file
 */

process INDEX {

  publishDir "${params.outdir}", mode: 'copy'
  input:
  path ref_genome

  output:
  path "index.*"

  script:
  """
  bowtie2-build $ref_genome index
  """
}

workflow {

    index_ch = INDEX(params.ref_genome)
}





