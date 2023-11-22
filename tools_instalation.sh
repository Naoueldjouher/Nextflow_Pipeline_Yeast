#!/bin/bash

# Install Nextflow
echo "Installing Nextflow..."
curl -s https://get.nextflow.io | bash
export PATH=$PATH:$(pwd)
nextflow -v

# Install Picard
echo "Installing Picard..."
PICARD_VERSION="2.25.6"  # Replace with the desired version
PICARD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"
mkdir -p tools/picard
curl -L -o tools/picard/picard.jar $PICARD_URL
echo "Picard installation complete."

# Install GATK
echo "Installing GATK..."
GATK_VERSION="4.2.0.0"  # Replace with the desired version
GATK_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip"
mkdir -p tools/gatk
curl -L -o tools/gatk/gatk.zip $GATK_URL
unzip tools/gatk/gatk.zip -d tools/gatk
chmod +x tools/gatk/gatk
echo "GATK installation complete."

# Install Trimmomatic
echo "Installing Trimmomatic..."
TRIMMOMATIC_VERSION="0.40"  # Replace with the desired version
TRIMMOMATIC_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
mkdir -p tools/trimmomatic
curl -L -o tools/trimmomatic/trimmomatic.zip $TRIMMOMATIC_URL
unzip tools/trimmomatic/trimmomatic.zip -d tools/trimmomatic
chmod +x tools/trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic.jar
echo "Trimmomatic installation complete."


# Install Samtools
echo "Installing Samtools..."
SAMTOOLS_VERSION="1.12"  # Replace with the desired version
SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
mkdir -p tools/samtools
curl -L -o tools/samtools/samtools.tar.bz2 $SAMTOOLS_URL
tar -jxf tools/samtools/samtools.tar.bz2 -C tools/samtools --strip-components 1
cd tools/samtools && make
export PATH=$PATH:$(pwd)
cd -
echo "Samtools installation complete."

# Install Bowtie2
echo "Installing Bowtie2..."
BOWTIE2_VERSION="2.4.4"  # Replace with the desired version
BOWTIE2_URL="https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"
mkdir -p tools/bowtie2
curl -L -o tools/bowtie2/bowtie2.zip $BOWTIE2_URL
unzip tools/bowtie2/bowtie2.zip -d tools/bowtie2
export PATH=$PATH:$(pwd)/tools/bowtie2
echo "Bowtie2 installation complete."

echo "All tools have been successfully installed."
