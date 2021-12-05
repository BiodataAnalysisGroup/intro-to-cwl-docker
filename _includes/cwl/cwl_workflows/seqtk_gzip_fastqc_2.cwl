#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}

inputs:
  # seqtk input
  fastqFile_array:
    type: File[]
  # # gzip file name output
  # seqtk_gzip_outName: 
  #   type: string

steps:
  subworkflow:
    run:
      class: Workflow
      
      inputs:
        fastqFile: 
          type: File
      
      outputs:
        # seqtk output
        seqtk_sample_output:
          type: File
          outputSource: seqtk_random_selection/seqtk_sample_out
        # seqtk gzip output
        seqtk_gzip_output:
          type: File
          outputSource: seqtk_gzip_step/gzip_out
        # output of fastqc step:
        zipFile:
          type: File
          outputSource: fastqc_quality_check/zipped_file
        htmlFile:
          type: File
          outputSource: fastqc_quality_check/html_file
        summaryFile:
          type: File
          outputSource: fastqc_quality_check/summary_file

      steps:
        # seqtk random subset selection (1%)
        seqtk_random_selection:
          run: ../cwl_tools/seqtk.cwl
          in: 
            fastq_file: fastqFile
          out: [seqtk_sample_out]
        # gzip step for seqtk output
        seqtk_gzip_step:
          run: ../cwl_tools/gzip.cwl
          in: 
            # outputFileName: seqtk_gzip_outName
            inputFile: seqtk_random_selection/seqtk_sample_out
          out: [gzip_out]
        # fastqc step
        fastqc_quality_check:
          run: ../cwl_tools/bio-cwl-tools/fastqc_2.cwl
          in: 
            reads_file: seqtk_gzip_step/gzip_out
          out: [zipped_file, html_file, summary_file]

    scatter: fastqFile
    in:
      fastqFile: fastqFile_array
    out: 
      - seqtk_sample_output
      - seqtk_gzip_output
      - zipFile
      - htmlFile
      - summaryFile

outputs:
  # seqtk output
  seqtk_sample_output2:
    type: 
      type: array
      items: File
    outputSource: subworkflow/seqtk_sample_output
  # seqtk gzip output
  seqtk_gzip_output2:
    type: 
      type: array
      items: File
    outputSource: subworkflow/seqtk_gzip_output
  # output of fastqc step:
  zipFile2:
    type: 
      type: array
      items: File
    outputSource: subworkflow/zipFile
  htmlFile2:
    type: 
      type: array
      items: File
    outputSource: subworkflow/htmlFile
  summaryFile2:
    type: 
      type: array
      items: File
    outputSource: subworkflow/summaryFile
