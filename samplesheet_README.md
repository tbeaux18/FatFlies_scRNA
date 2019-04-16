# FatFlies_scRNA
Single cell RNA-seq Analysis Pipeline
Sample Sheet Tutorial

Important: all files need to need be in the same folder

Date: Input the date that you are running test data
Run_name : This can be something like 
Library_Prep: this should remain cel-seq2, unless you begin to use a new library prep method
Basename: choose a name that all files will be created with for this pipeline
Fastq_Read1: put the absolute path to the read1 fastq file
Fastq_Read2: put the absolute path to the read2 fastq file
Ref_genome: put the absolute path to the .fa file 
Annotation: put the absolute path to the .gtf file

[ZUMI}
Bc_filter_num_bases 1			
bc_filter_phred	20			
bc_ham_dist	1			
umi_filter_num_bases	1			
umi_filter_phred	20			
umi_ham_dist	1			
zum_start_stage	Filtering

These values are taken from the zUMI page and should be treated like defaults

[DIFF_EXP]

test_group
13
If there are 4 groups, you write the testing groups as the integer. You do not need to write a comma or space between the group numbers. 

control_group
24
If there are 4 groups, you write the testing groups as the integer. You do not need to write a comma or space between the group numbers. 

 
[ADAPTERS]				
Adapter	TGGAATTCTCGG	GATCGTCGGACT		

Put the 3' and 5' adapters in the cells

[DATA]				
cell_id	sample_name	barcode_index	barcode_sequence	experiment_group
 1	        F1.2	        cel-4        	AGCTTC	            1

cell_id numbers each sample that is put into the sample sheet, starting at 1 and growing in incraments of 1
sample_name is F for fed and S for starved and labels 1 for group number and the 2 for cell from the grouping. ie Fed group 1, cell 2
barcode_seqeunce is the sequence that is returned from the library prep
experiment_group is in reference to the experiment. 1,2 are for the fed and 3,4 are for starved. 
