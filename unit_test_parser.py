#!/usr/bin/env python
# coding: utf-8

# In[1]:
#4/9/19
#Amanda 

#BN base name, RN run name 
def unit_base_name(RN,BN):
    assert RN.startswith(BN), "Your run name and base name do not share the same base."

def unit_library_prep(LP):   
    assert LP == "cel-seq2" , "you should have typed cel-seq2"
      
def unit_bc_filter_num_bases(NB):
    assert NB == 1, "bc_filter_num_bases should have been 1."
    
def unit_bc_filter_phred(BFP):
    assert BFP == 20 , "bc_filter_phred should be 20"

def unit_bc_ham_dist(HD):
    assert HD == 1 , "bc_ham_dist should be 1"

def unit_umi_filter_num_bases(UNM):
    assert UNM == 1 , "umi_filter_num_bases should be 1"

def unit_umi_filter_phred(UFP):
    assert UFP == 20, "umi_filter_phred should be 20"

def unit_umi_ham_dist(UHD):
    assert UHD == 1, "umi_ham_dist should be 1"
#making sure stage is the right part
def unit_zum_start_stage(ZSS):
    assert ZSS == "Filtering" , "zum_start_stage should be 'Filtering'"

#can use for 3 and 5
#also i am not 100% sure this is a good test
def unit_adaptor(UA):
    assert len(UA)== 12, "Your adaptor should be 6 base pairs long"

#testing for read 1 and read 2
#making sure files are the right kind
def unit_fastq(UF):
    assert UF.endswith("fastq"), "this needs to be a .fastq file"

#making sure files are the right kind
def unit_ref_genome(RF):
    assert RF.endswith("fa"), "this needs to be a .fa file"
#making sure files are the right kind
def unit_annotation(UA):
    assert UA.endswith("gtf"), "this needs to be a .gtf file"

#def unit_test_groups(TG):
    
#ALL TESTS UNDER HERE
#run name goes first and base name goes second to make sure things are all matched up
#unit_base_name("this_test" , "this")
#unit_library_prep('cel-seq2')
#unitbc_filter_num_bases(1)
#unit_bc_filter_phred(20)
#unit_bc_ham_dist(1)
#unit_umi_filter_num_bases(1)
#unit_umi_filter_phred(20)
#unit_umi_ham_dist(1)
#unit_zum_start_stage('Filtering')
#3 and 5 here 
#unit_adaptor("GGGGGGGGGGGG")
#read 1 and 2
#unit_fastq("this.fastq")
#unit_ref_genome("this.fa")
#unit_annotation("this.gtf")
