
</body>
</html>

<h1>Integrative Analysis Pipeline</h1>

<p>This is a Nextflow pipeline for RNA-Seq analysis, covering the steps of SRA downloading, FASTQ trimming, alignment, Picard indexing, and GATK variant calling.</p>

<h2>Prerequisites</h2>

<ul>
  <li><a href="https://www.nextflow.io/">Nextflow</a></li>
  <li>Other dependencies such as Trimmomatic, BWA, Picard, GATK</li>
</ul>

<h2>Usage</h2>

<p>Clone the repository:</p>
<code>git clone https://github.com/your-username/your-repo.git</code>

<p>Navigate to the project directory:</p>
<code>cd your-repo</code>

<p>Run the pipeline:</p>
<code>nextflow your_pipeline.nf --srr_numbers </code>
<h2>Workflow</h2>

<p>The pipeline consists of the following processes:</p>

<ul>
  <li><strong>downloadsrr.nf:</strong> Downloads SRA data using SRA Toolkit</li>
  <li><strong>trimming.nf:</strong> Trims FASTQ files using Trimmomatic</li>
  <li><strong>alignment.nf:</strong> Aligns trimmed reads using BWA (example)</li>
  <li><strong>sortandpicard.nf:</strong> Indexes aligned BAM files using Picard</li>
  <li><strong>Gatk.nf:</strong> Calls variants using GATK</li>
</ul>

<h2>Parameters</h2>

<p>Customize the pipeline by modifying the parameters in the script:</p>

<!-- Step 1: DOWNLOAD -->
<h2>Step 1: SRA Data Download (DOWNLOAD)</h2>
<p>Downloads raw sequencing data from SRA using the SRA Toolkit.</p>
<p><strong>Parameters:</strong></p>
<code>params.srr_numbers</code> - SRA accession number to download.
<code>params.outdir</code> - Output directory for downloaded files.
<code>params.sratoolkit</code> - Path to the SRA Toolkit executable.

<!-- Step 2: Trimming (runTrimmomatic) -->
<h2>Step 2: Quality Trimming (runTrimmomatic)</h2>
<p>Trims raw FASTQ files for quality using Trimmomatic.</p>
<p><strong>Parameters:</strong></p>
<code>params.trim_path</code> - Path to Trimmomatic installation.
<code>params.trim_jar</code> - Path to Trimmomatic JAR file.
<code>params.java_path</code> - Path to Java executable.
<code>params.fastqc</code> - Path to FastQC executable.
<code>params.srr_numbers</code> - SRA accession number for trimming.
<code>params.outdirtrimm</code> - Output directory for trimmed files.
<!-- Step 3: Bowtie2 Alignment (bowtieAlignment) -->
<h2>Step 3: Bowtie2 Alignment</h2>
<p>Aligns trimmed and merged FASTQ files to a reference genome using Bowtie2.</p>
<p><strong>Parameters:</strong></p>
<code>params.outdiralign</code> - Output directory for alignment results.
<code>params.bowtie2_exec</code> - Path to Bowtie2 executable.
<code>params.index_path</code> - Path to the Bowtie2 index for the reference genome.
<code>params.samtools_path</code> - Path to Samtools executable.
<code>params.threads</code> - Number of threads for parallel processing.
<code>params.srr_numbers</code> - SRA accession number for alignment.
<code>params.outdirtrimm</code> - Output directory for trimmed and merged files.

<!-- Bowtie2 Alignment Process -->
<code>process bowtieAlignment</code>
<pre>
  <code>
  ${params.bowtie2_exec} \\
      -x "${params.index_path}index" \\
      -U "${params.outdirtrimm}" \\
      -S "${srr_numbers}.sam"
  </code>
</pre>

<!-- SAM to BAM Conversion (samtobam) -->
<h2>Step 4: Convert SAM to BAM</h2>
<p>Converts SAM files to BAM format using Samtools.</p>
<p><strong>Parameters:</strong></p>
<code>params.outdiralign</code> - Output directory for alignment results.

<!-- SAM to BAM Conversion Process -->
<code>process samtobam</code>
<pre>
  <code>
  ${params.samtools_path} view -bS $sam > "output.bam"
  </code>
</pre>
<!-- Step 5: Picard/GATK Processing -->
<h2>Step 5: Picard/GATK Processing</h2>
<p>Utilizes Picard and GATK tools for additional processing of aligned BAM files.</p>
<p><strong>Parameters:</strong></p>
<code>params.samtools_path</code> - Path to Samtools executable.
<code>params.threads</code> - Number of threads for parallel processing.
<code>params.srr_numbers</code> - SRA accession number for processing.
<code>params.align</code> - Output directory for alignment results.
<code>params.picard</code> - Path to Picard JAR file.
<code>params.javaCmd</code> - Java command for running Picard and GATK.
<code>params.ref</code> - Path to the reference genome FASTA file.
<code>params.fai</code> - Path to the reference genome FASTA index file.
<code>params.tbi</code> - Path to the reference genome VCF index file.
<code>params.dict</code> - Path to the reference genome dictionary file.
<code>params.vcf</code> - Path to the reference genome VCF file.
<code>params.picard_gatk</code> - Output directory for Picard/GATK results.

<!-- Sort BAM files (sortBam) -->
<code>process sortBam</code>
<pre>
  <code>
  ${params.samtools_path} sort -o ${params.srr_numbers}_output.sorted.bam ${reads}
  </code>
</pre>

<!-- Add Read Groups (addReadGroupProcess) -->
<code>process addReadGroupProcess</code>
<pre>
  <code>
  ${params.javaCmd} ${params.picard} AddOrReplaceReadGroups \
    I=${bamfile} \
    O=${params.srr_numbers}_output.sorted.RD.RG.bam \
    RGID=1 RGLB=${params.srr_numbers} RGPL=ILLUMINA RGPU=run RGSM=${params.srr_numbers}
  </code>
</pre>






</body>
</html>
