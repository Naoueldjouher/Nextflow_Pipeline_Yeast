params.url = "ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/"
params.base = "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome."
params.chrs = (1..16).collect { it.toString() }
params.downloadDir = "downloaded_chromosome_files"
params.outputFasta = "Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz"

// Define the process for downloading a chromosome sequence
process downloadChr {
    input:
    val chr

    output:
    file("${params.downloadDir}/${params.base}${chr}.fa.gz")

    script:
    """
    mkdir -p ${params.downloadDir}
    curl -o ${params.downloadDir}/${params.base}${chr}.fa.gz ${params.url}${params.base}${chr}.fa.gz
    """
}

// Define the process for concatenating and compressing the FASTA files
process concatenateAndCompress {
    input:
    file downloadedChromosomes

    output:
    file params.outputFasta

    script:
    """
    cat ${downloadedChromosomes} > ${params.outputFasta}
    bgzip -k ${params.outputFasta}
    """
}

// Define the process for indexing the reference genome with samtools faidx
process indexReference {
    input:
    file params.outputFasta

    script:
    """
    samtools faidx ${params.outputFasta}
    """
}

// Define the process for creating a sequence dictionary with Picard
process createSequenceDictionary {
    input:
    file params.outputFasta

    output:
    file("${params.outputFasta}.dict")

    script:
    """
    picard CreateSequenceDictionary \
        R=${params.outputFasta} \
        O=${params.outputFasta}.dict
    """
}

// Define the process for creating a Tabix index with GATK
process createTabixIndex {
    input:
    file params.outputFasta

    script:
    """
    gatk IndexFeatureFile \
        -F ${params.outputFasta}
    """
}

workflow {
    // Define the download process for each chromosome in parallel
    downloadedFiles = downloadChr(params.chrs)

    // Define the concatenation and compression process after all downloads are completed
    final_genome = concatenateAndCompress(downloadedFiles)

    // Index the reference genome with samtools faidx
    indexReference(final_genome)

    // Create a sequence dictionary with Picard
    createSequenceDictionary(final_genome)

    // Create a Tabix index with GATK
    createTabixIndex(final_genome)
}
