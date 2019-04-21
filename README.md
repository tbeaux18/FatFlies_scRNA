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

On the local machine:
  * Create a directory that contains the following:
      * Reference_Genome_FASTA file
      * Reference_Annotation gtf file
      * Compressed Read 1 FASTQ (gz)
      * Compressed Read 2 FASTQ (gz)
      * SampleSheet.csv
  * Refer to this repository's wiki on creating the sample sheet
  * Sample sheet **must** be named SampleSheet.csv

After creating directory and adding the files, change into it and clone this repository.
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

Running this command will trigger the docker images to build, and once built, the pipeline will run a docker container which will run the pipeline. The docker images were built to run mostly on Linux/MacOS, and theoretically run on Windows. This pipeline could break if put on AWS due to its own dependencies. 


### Prerequisites

Refer to the Sample Sheet Wiki for creating the Sample Sheet
The only required software is Docker, which the latest installation documentation can be found at this link:
  * https://docs.docker.com/v17.12/install/



### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

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

