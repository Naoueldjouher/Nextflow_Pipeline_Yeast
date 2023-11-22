params {
    url = "ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/"
    base = "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome."
    chrs = (1..16).collect { it.toString() }
    downloadDir = "downloaded_chromosome_files"
    outputFasta = "Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz"
}

// Define the process for downloading a chromosome sequence
process downloadChr {
    input:
    val chr

    output:
    file("${downloadDir}/${base}${chr}.fa.gz") into downloadedFiles

    script:
    """
    mkdir -p ${downloadDir}
    curl -o ${downloadDir}/${base}${chr}.fa.gz ${params.url}${params.base}${chr}.fa.gz
    """
}

// Define the process for concatenating and compressing the FASTA files
process concatenateAndCompress {
    input:
    file downloadedChromosomes from downloadedFiles

    output:
    file params.outputFasta

    script:
    """
    cat ${downloadedChromosomes} > ${params.outputFasta}
    bgzip -k ${params.outputFasta}
    """
}

// Run the download process for each chromosome in parallel
downloadChr(params.chrs)

// Run the concatenation and compression process after all downloads are completed
concatenateAndCompress(downloadedFiles)
