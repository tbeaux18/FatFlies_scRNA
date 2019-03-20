# CELSEQ2 Workflow using STAR aligner
This portion of the pipeline can be carried out in four steps
1. Build STAR index (done)
2. Make config.yaml file (done)
3. Create an experiment table (done)
4. Run the celseq2 pipeline (in progress)

## Commands for building the STAR index (complete)
    mkdir dm6_ERCC_GAL4_GFP_STAR

    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir dm6_ERCC_GAL4_GFP_STAR \
    --sjdbGTFfile dm6.refGene.CEL-Seq.gtf --sjdbOverhang 49 --genomeFastaFiles dm6_ERCC_GAL4_GFP.fa

    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir dm6_ERCC_GAL4_GFP_STAR_noGTF \
    --sjdbOverhang 49 --genomeFastaFiles dm6_ERCC_GAL4_GFP.fa

## Commands for running CELSEQ2 (work in progress)
First, the configuration and experiment table must be initialized.

    new-configuration-file -o config.yaml
    new-experiment-table -o experiment_table.txt

These files were edited to contain relevant parameters and file paths. See config.yaml and experiment_table.txt for more details. Finally, the celseq2 command is run.

    celseq2 --config-file /home1/jjackson10/FatFlies_scRNA/celseq2/config.yaml \
    --experiment-table /home1/jjackson10/FatFlies_scRNA/celseq2/experiment_table.txt \
    --output-dir /home1/jjackson10/FatFlies_scRNA/celseq2/results \
    --dry-run \
    -j 16

## Running away from server
    screen
    ...
    Ctrl-a Ctrl-d
    screen -list
    screen -r
