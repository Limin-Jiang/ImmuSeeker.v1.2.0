# ImmuSeeker1.2.0
This project encompasses all the code utilized in the paper titled "Deep mining of immune-related gene family signatures through lineage reconstruction"

## Paper Description
Immune gene family signatures, such as HLA and Killer-cell Immunoglobulin-like Receptor (KIR), can be extracted from sequencing data. Previous methods offer limited scope, functionality, and preservation. By examining the phylogenetic structures of HLA and KIR, we developed ImmuSeeker for accurate HLA/KIR gene and allele and protein sequence identification, assessment of HLA/KIR gene allele-specific, and protein-sequence-specific expression, diversity, Bayesian-based zygosity inference, and contrastive neural network-based HLA/KIR haplotype group comparison. 

![Framework](https://github.com/Limin-Jiang/HLA_autoimmune/blob/main/Figure.JPG)


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example](#Example)
- [Additional Information](#Additional)
- [Contact](#contact)

## Installation ImmuSeeker

### Via docker (Recommended)

```bash
docker pull lxj423/immuseeker1.2.0
```

### Manual installation

####  Prerequisites

- Samtools (Copyright (C) 2008-2024 Genome Research Ltd.)
```bash
apt -y install samtools
```

- bowtie

```bash
apt -y install bowtie
```

- R  >=v3.6 ( With packages: data.table, dplyr, ggraph, igraph, docopt, phyloseq, stringr)
```bash
apt -y install r-base
Rscript -e "install.packages(c('data.table', 'dplyr', 'ggraph', 'igraph', 'docopt','stringr','Biostrings','phyloseq'), repos='https://cran.r-project.org')"
```

####  Install
To get started with this project, follow these steps:

```bash
# Clone the repository
git clone https://github.com/Limin-Jiang/ImmuSeeker.v1.2.0.git

# Navigate to the project directory
cd ImmuSeeker
export PATH=/your/directory/ImmuSeeker.v1.2.0:$PATH
```



## Usage

```bash
Usage: ImmuSeeker [-HLA | -KIR] [-SwithB | -SwithF | -Update] -i <Parameter1> -v <Parameter2> -c <Parameter3> -n <Parameter4>  -p <Parameter5> -p1 <Parameter6> -o <Parameter7>  -er <Parameter8> -ex <Parameter9> -pt <Parameter10> -dv <Parameter11> -nr <Parameter12> -ca <Parameter13> -cov <Parameter14> -r <Parameter15> -q <Parameter16> -rv <Parameter17>
Options:
    Required:
	-HLA     Invoke the HLA calling process.
	-KIR     Invoke the KIR calling process.
	-SwithB  Input alignment BAM file.
	-SwithF  Input FASTQ file.
	-Update  Update to latest IMGT-HLA/KIR version.
	-i       Set the input directory and filename. Example: 'your/input/directory/inputfile.bam'.
		 If using -SwithF, specify FASTQ file(s) with '-i file1,file2' (two files) or '-i file' (one file).
    Optional:
	-v       If using -SwithB, specify the genome version with '-v hg37' or '-v hg38'.
	-c       Specify min number of supported unique reads (default: 1).
	-ca      Specifies the minimum number of supporting reads required for each allele to be considered (default: 100).
	-cov     Specifies the minimum coverage threshold required for an allele to be retained (default: 98%).
	-r       Specifies the minimum ratio threshold required between two alleles for both to be retained (default: 50%).
	-q       Specifies the minimum quality threshold for reads to be included in the analysis (default: 30).
	-n       Max mismatches in alignment. (default: 0).
	-p       Set a noninformative flat prior to allow the data to have a strong influence on the posterior distribution. (default: -p '(1/3,1/3,1/3)').
	-p1      Set the probability of observing an allele in genotype (default: 1/2).
	-o       Set the output directory and filename. Example: 'your/output/directory/outfile'.
	-ex      Specify whether to output gene expression values(default: -ex false).
	-pt      Specify whether to output phylogenetic tree for HLAs (default: false).
	-dv      Specify whether to output diversity analysis results based on the number of unique reads ('dvr') or gene expression ('dve'). (default:  false).
	-nr      Specify the type of reads to account for when calculating gene expression. Choose between HLA reads ('HLAn') or total reads ('Totaln'). Calculating using total reads requires more time. (default: 'HLAn').
	-er      Specify the sequencing error ratio (default: 0.02).
	-rv      Specify the input HLA/KIR reference file name. (Default: uses the latest available version).
	--help   Display this help message.
```
## Example

#### Example based on docker

cd /your/data/folder/ (that directory includes your data folder/file)
##### Example for updating to the latest IMGT-HLA/KIR version.
```bash
docker run -v /your/data/folder:/ImmuSeeker_data -it immuseeker1.2.0 -KIR -Update 
docker run -v /your/data/folder:/ImmuSeeker_data -it immuseeker1.2.0 -HLA -Update
```
If you prefer to install a specific version of the reference, you must first download the hla_nuc.fasta or kir_nuc.fasta files from the IPD-IMGT/HLA Database. Once downloaded, run the following commands to install the reference files:
```bash
sudo docker run -v /your/data/folder:/ImmuSeeker_data -it immuseeker1.2.0 -KIR -Update -rv ../ImmuSeeker_data/kirla_nuc.fasta
sudo docker run -v /your/data/folder:/ImmuSeeker_data -it immuseeker1.2.0 -HLA -Update -rv ../ImmuSeeker_data/hla_nuc.fasta
```
The most recent HLA/KIR reference files will be stored in /your/data/folder for future use.

##### Next, ensure that your BAM or FASTQ files are stored in /your/data/folder. For example, if Example.bam is located in /your/data/folder, execute the following command:
```bash
sudo docker run -v /your/data/folder:/ImmuSeeker_data -it immuseeker1.2.0 -HLA -SwithB -i bamfile.bam -ca 40  -cov 95 -ex true -pt true -dv dve -q 30
```


#### Example in local

##### Example for updating to the latest IMGT-HLA/KIR version.
```bash
./ImmuSeeker -HLA -Update 
./ImmuSeeker -KIR -Update
```
If you prefer to install a specific version of the reference, you must first download the HLA_nuc.fasta or kir_nuc.fasta files from the IPD-IMGT/HLA Database. Once downloaded, run the following commands to install the reference files:

```bash
./ImmuSeeker -HLA -Update -rv '/your/data/folder/HLA_nuc.fasta'
./ImmuSeeker -KIR -Update -rv '/your/data/folder/kir_nuc.fasta'
```
The installed HLA/KIR reference files will be stored in /ImmuSeeker/data, replacing the previous reference files for future use.

##### This is an example for input a bam file to detect HLA:
```bash
./ImmuSeeker -HLA -SwithB -i Example.bam -o Example_output -n 0 -c 0 -v hg38 -cov 90 -q 20
```

##### This is an example for input one or two fq file to detect HLA:
```bash
./ImmuSeeker -HLA -SwithF -i Example.fq -er 0.02 -n 0 -c 0 -o output -ex true -dv dve -pt true
```

or

```bash
./ImmuSeeker -HLA -SwithF -i file1.fq,file2.fq -er 0.02 -n 0 -c 0 -o output -ex true
```

This is an example for detecting KIR:
```bash
./ImmuSeeker -KIR -SwithB -i Example.bam -o Example_output
./ImmuSeeker -KIR -SwithF -i Example.fq -o Example_output
```


## Output File Descriptions

Upon completion, several output files will be generated, each providing specific information related to HLA/KIR allele detection and analysis:

### output.all.csv: 

Contains all potential allele candidates. This file includes detected allele IDs, gene names, the number of unique supporting reads, and associated allele levels. If the -ex true option is enabled, it will also include allele expression values.

### output.candidate.csv: 

Contains the final set of predicted alleles. Similar to output.all.csv, it includes allele IDs, gene names, unique supporting read counts, and allele levels. Expression values are included if -ex true is specified.

### output_genetype.csv: 

Reports the inferred zygosity of each gene, estimated using Bayesian modeling.


### output_diversity_ex.csv/output_diversity_re.csv: 

Provide results from diversity analysis, either based on gene expression (output_diversity_ex.csv) or read counts (output_diversity_re.csv).

### output_evolution_graphs.pdf: 

Displays phylogenetic trees representing the evolutionary relationships among the detected HLA alleles.


## Additional information

This project also includes an analysis of HLA diversity, as well as a comparison of HLA haplotype groups using a Contrastive Neural Network approach. The corresponding code and example datasets are available in the Example_Analysis_CL_and_diversity directory within the original pipeline, accessible at <a href="https://github.com/Limin-Jiang/ImmuSeeker">ImmuSeeker </a>.


## Contact

Any Comments or questions, Please contact:

Yan Guo, Ph.D, yanguo1978@gmail.com

Limin Jiang, Ph.D, lxj423@miami.com
