#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
# input of trim_galore step:  
  raw_fastq_file_1: 
    type: File
  raw_fastq_file_2: 
    type: File
# input for bwa mem
  ref_genome_2nd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
# input for samtools view sam2bam
  sort_type:
    type: boolean
    default: false
# input for bcftools mpileup
  ref_genome_3rd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .fai
  output_file_type:
    type: string
    default: "u"
# input for bcftools call
  output_file_type_2nd:
    type: string
    default: "v"
  type_of_ploidy: 
    type: int
  output_file_name:
    type: string
  type_of_model:
    type: string
    default: "m"

outputs: 
# output of trim_galore step:
  fastq_1_trimmed_file:
    type: File
    outputSource: trimming_reads/fastq1_trimmed
  fastq_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed
  fastq_unpaired_1_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq1_trimmed_unpaired
  fastq_unpaired_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed_unpaired
  trim_galoreLog:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trim_galore_log
  trimmed_fastqc_html_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_html
  trimmed_fastqc_zip_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_zip
# output for bwa mem
  bwa_mem_file:
    type: File
    outputSource: bwa_mem_mapping_reads/reads_stdout
# output for samtools view sam2bam
  sam_to_bam_file:
    type: File
    outputSource: sam_to_bam_conversion/bam
# output for samtools sort
  bam_sorted_file:
    type: File
    outputSource: bam_file_sorting/bam_sorted
# output for bcftools mpileup call pipeline
  bcftools_vcf_file:
    type: File
    outputSource: bcftools_call_step/output_vcf_file


steps:
# trim_galore trimming of reads
  trimming_reads:
    run: ../cwl_tools/bio-cwl-tools/trim_galore.cwl
    in:
      fastq1: raw_fastq_file_1
      fastq2: raw_fastq_file_2
    out:
      - fastq1_trimmed 
      - fastq2_trimmed
      - fastq1_trimmed_unpaired 
      - fastq2_trimmed_unpaired 
      - trim_galore_log
      - trimmed_fastqc_html
      - trimmed_fastqc_zip
# bwa mem mapping reads:
  bwa_mem_mapping_reads:
    run: ../cwl_tools/bio-cwl-tools/BWA-Mem.cwl
    in: 
      InputFile_1: trimming_reads/fastq1_trimmed # trimmed_reads_1
      InputFile_2: trimming_reads/fastq2_trimmed
      Index: ref_genome_2nd
    out: [reads_stdout]
# SAM_to_BAM_conversion:
  sam_to_bam_conversion:
    run: ../cwl_tools/bio-cwl-tools/samtools_view_sam2bam.cwl
    in:
      sam: bwa_mem_mapping_reads/reads_stdout # sam_file
    out: [bam]
# BAM_sorting:
  bam_file_sorting:
    run: ../cwl_tools/bio-cwl-tools/samtools_sort.cwl
    in:
      bam_unsorted: sam_to_bam_conversion/bam  #bam_file
      by_name: sort_type
    out: [bam_sorted]
# bcftools_mpileup_step:
  bcftools_mpileup_step:
    run: ../cwl_tools/bcftools_mpileup.cwl
    in:
      refGenome: ref_genome_3rd
      outputFileType: output_file_type
      inputBam: bam_file_sorting/bam_sorted # bam_file_2nd
    out: [pileupOut]
# bcftools_call_step:
  bcftools_call_step:
    run: ../cwl_tools/bcftools_call.cwl
    in:
      outputFileType: output_file_type_2nd
      ploidyType: type_of_ploidy
      outputFileName: output_file_name
      inputFileMpileup: bcftools_mpileup_step/pileupOut
      calling_model: type_of_model
    out: [output_vcf_file]

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
