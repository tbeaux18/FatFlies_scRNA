# FatFlies_scRNA

### Single cell RNA-seq Analysis Pipeline

This pipeline is a wrapper of various tools associated with single cell RNA-seq analysis. The pipeline is optimized for 3' based library preparation methods in NGS for single cell. This pipeline is specifically optimized for CEL-Seq2 libraries that have pseudo-paired reads. I plan to generalize this pipeline for expanded use.

Current library prep requirements are:
  * CEL-Seq2
  * Already lane merged fastq files that result in a single pair of fastqs.
    * Read 1 FASTQ
      * The UMI in position 1 - 6
      * The cell barcode in position 7-12
    * Read 2 FASTQ
      * The transcript 1-50 (or the max length of the second record)
  * Sequenced on Illumina platform
  * Relatively low throughput due to sample sheet requirements
      

## Getting Started

The pipeline is dockerized and contains all the necessary software dependencies to run quality control on raw reads, perform alignment and UMI counting, and run a basic differential expression analysis with a simple experiment design.

### Prerequisites

Refer to the Sample Sheet Wiki for creating the Sample Sheet
The only required software is Docker, which the latest installation documentation can be found at this link:
  * https://docs.docker.com/v17.12/install/
  
### Host machine directory setup

On the local machine:
  * Create a directory that contains the following:
      * Reference_Genome_FASTA file
      * Reference_Annotation gtf file
      * Compressed Read 1 FASTQ (gz)
      * Compressed Read 2 FASTQ (gz)
      * SampleSheet.csv
  * Refer to this repository's wiki **Sample Sheet README** on creating the sample sheet
  * Sample sheet **must** be named SampleSheet.csv

After creating directory and adding the files, change into the directory and clone this repository.
```
cd pipeline_example
git clone https://github.com/tbeaux18/FatFlies_scRNA.git
```
This results in the following directory structure.
```
root:pipeline_example example$
.
├── FatFlies_scRNA
├── SampleSheet.csv
├── read1.fastq.gz
├── read2.fastq.gz
├── ref_genome.fa
└── ref_genome_annotation.gtf
```

To run the pipeline:
```
cd FatFlies_scRNA
./run_pipeline.sh
```

### Pipeline Procedure
By running this command:
  * Docker images will build
  * Docker runs the **ubuntur35:pipeline** container
  * Docker mounts the parent directory of FatFlies_scRNA
  * Sample Sheet is parsed and exchanges information with necessary input files such as zUMIs.yaml file
  * Quality Control of Raw Data
  * STAR Index is built with the provided FASTA, the GTF is not included in the index build
  * zUMIs pipeline begins
    * Stages
      * Filtering
      * Mapping
      * Counting
    * Output directory
      * <sample_sheet_basename>_zumi_output
        * BAM Files
        * zUMI log files
        * Count matrices
        * Alignment plots
  * Differential Expression Rscript runs
    * Performs DE on respective Control vs. Test Groups
    * Performs DE on respective Condition1 vs Condition2 Test Groups 
    * Outputs **r_results** Directory
      * HTML files for interactive scatterplots for each comparison
      * CSV files for all and significant differential expresssed genes for each comparison
  * Logs are outputted for each step and are located in ~/Fatflies_scRNA/logs directory



## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Authors

* **Timothy Baker**
* **Amanda Wills**
* **Jasen Jackson**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Dr. Jennifer Beshel, PhD
* Dr. Catherine Putonti, PhD
* Loyola University Chicago
* CEL-Seq2 authors

