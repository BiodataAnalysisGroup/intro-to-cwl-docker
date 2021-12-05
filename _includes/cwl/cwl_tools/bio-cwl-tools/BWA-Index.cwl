#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"

baseCommand: [bwa, index]

inputs:
  InputFile:
    type: File
    format: edam:format_1929  # FASTA
    inputBinding:
      position: 200
    
  IndexName:
    type: string
    inputBinding:
      prefix: "-p"
      position: 1
      # valueFrom: $(self).bwt

#Optional arguments

  algoType:
    type: 
      - "null"
      - type: enum
        symbols:
        - is
        - bwtsw
    inputBinding:
      prefix: "-a"

outputs:
  index:
    doc: all secondary files required for bwa mem
    type:
      type: array
      items: File
    outputBinding:
      glob: "$(inputs.IndexName)*"

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
