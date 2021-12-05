#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.14--h88f3f91_0"

baseCommand: [bcftools, mpileup]

inputs:
  refGenome:
    type: File
    format: edam:format_1929 # FASTA
    inputBinding:
      prefix: "-f"
      position: 1
    # secondaryFiles: # .fai is created automatically when bcftools mpileup is called...
    #  - .fai 
  
  outputFileType: # select type of output file: b|u|z|v
    type: string
    inputBinding:
      position: 2
      prefix: -O
      separate: false
    default: "u"
  
  inputBam:
    type: File
    inputBinding:
      position: 100

outputs:
  pileupOut:
    type: stdout

# stdout: mpileup_output.txt

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
