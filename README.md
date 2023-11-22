
</body>
</html>

<h1>Integrative Analysis Pipeline</h1>
<h2>Project Introduction</h2>

<p>Welcome to the Nextflow Yeast Analysis Pipeline!</p>

<p>This project aims to reproduce a previously conducted study on Chip-seq analysis of the Yeast Genome. The original study can be found (https://github.com/Naoueldjouher/Chip-seq_Yeast_Genome/blob/main/README.md) for further information about the study design, methodologies, and results.</p>

<!-- Reproduction with Nextflow -->
<h2>Reproducing the Pipeline with Nextflow</h2>

<p>In this project,the goal is to learn how to recreate the analysis in a more scalable, modular, and reproducible manner using Nextflow.</p>

<p>By following the steps outlined in this README, you will be able to execute the Nextflow pipeline and compare the results with the original study. This serves as a valuable learning experience in workflow management and reproducible research.</p>



<p>This is a Nextflow pipeline for RNA-Seq analysis, covering the steps of SRA downloading, FASTQ trimming, alignment, Picard indexing, and GATK variant calling.</p>

<h2>Prerequisites</h2>

<ul>
  <li><a href="https://www.nextflow.io/">Nextflow</a></li>
  <li>Other dependencies such as Trimmomatic, BWA, Picard, GATK</li>
</ul>

<h2>Clone the Repository</h2>

<p>To get started with the Nextflow Yeast Analysis Pipeline, clone the repository using the following command:</p>
<code>git clone https://github.com/Naoueldjouher/Nextflow_Pipeline_Yeast</code>

<p>Navigate to the project directory:</p>
<code>cd Nextflow_Pipeline_Yeast</code>

<!-- Usage -->
<h2>Usage</h2>
<p>Once you have cloned the repository and navigated to the project directory, you can run the pipeline with the following command:</p>
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

<h2>The Size Discrepancy</h2>

<p>During the execution of this pipeline, a notable observation was made regarding the size of the output files after trimming. Specifically, for example the trimmed paired-end FASTQ files for  <code>SRR1811834</code> exhibit a substantial difference in size:</p>

<ul>
  <li><code>SRR1811834_R1_paired_trimmed.fastq.gz</code>: 605 MB</li>
  <li><code>SRR1811834_R2_paired_trimmed.fastq.gz</code>: 44 MB</li>
</ul>

<p>This discrepancy prompts further investigation into the factors contributing to the differing sizes between the forward (R1) and reverse (R2) reads after the trimming process.</p>
<body>
<!DOCTYPE html>
<html lang="en">



<body>
<h2>Trimming Discrepancy Investigation</h2>
 <p>While executing the trimming step of this RNA-Seq analysis pipeline, a notable difference in output file sizes was observed between Nextflow and standalone execution. This section explores hypotheses and provides suggested actions to investigate and address the observed discrepancies.</p>
    <h3>Next Steps:</h3>

  <ol>
        <li><strong>Detailed Logging:</strong> Enabled detailed logging in both Nextflow and standalone executions. Examined the log files to identify any specific differences in the execution environment, parameters, or tool versions.</li>
        <li><strong>Intermediate File Inspection:</strong> Investigated intermediate files generated during the trimming process in both scenarios. Compare the content and sizes of these intermediate files to identify any discrepancies.</li>
        <li><strong>Tool Benchmarking:</strong> Benchmark the performance of the trimming tool outside Nextflow under similar conditions as the Nextflow workflow. Compares the runtime, resource utilization, and output file sizes.</li>
        <li><strong>Workflow Isolation:</strong> Isolated the trimming step from the Nextflow workflow and executed it independently. Compared the results to ensure that Nextflow itself is not introducing unintended modifications.</li>
        <li><strong>Tool Parameter Tuning:</strong> Experiment with different parameters of the trimming tool, both within and outside Nextflow. Adjust parameters related to quality thresholds, adapter removal, and other relevant settings.</li>
    
  </ol>

  <p>Given the lack of conclusive information from the initial investigations, a decision was made to proceed with the downstream steps of the pipeline using only <code>Pair_1</code>.</p>

<h2>GATK ApplyBSQR Issue</h2>

<p>During the execution of this pipeline, an issue was encountered with the GATK ApplyBSQR step. Previously, the same sample was processed using this step without any errors. However, in the current run using Nextflow, the pipeline reported an error indicating a mismatch between the reads and the reference genome. The reference genome used in both runs was the same, and the input reads were identical.</p>

<p><strong>Error Message:</strong></p>
<pre>
Error: There was an error in ApplyBSQR. No mismatches found between the read and the reference genome.
</pre>

<p><strong>Investigation Steps:</strong></p>

<ul>
  <li>Checked the reference genome: Confirmed that the reference genome used in the pipeline matched the one used successfully in previous runs.</li>
  <li>Validated the reads: Used Picard's ValidateSamFile and Samtools' quickcheck to ensure the integrity of the input reads. No errors were reported.</li>
  <li>Nextflow Execution: Reviewed Nextflow execution logs and parameters to identify any potential issues in the workflow.</li>
</ul>

<p><strong>Resolution:</strong></p>

<p>Despite the thorough checks, the exact cause of the ApplyBSQR error remains unclear. To address this issue and proceed with the analysis, a temporary solution was implemented by bypassing the ApplyBSQR step in the current run. This decision was made based on the assurance that the input data (reads and reference genome) were consistent and valid.</p>

<p><strong>Next Steps:</strong></p>

<p>The ApplyBSQR step will be revisited in subsequent pipeline updates, and further investigation will be conducted to identify the root cause of the error.</p>

</body>

</html>








