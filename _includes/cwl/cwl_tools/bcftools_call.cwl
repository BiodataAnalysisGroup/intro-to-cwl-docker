#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.14--h88f3f91_0"

baseCommand: [bcftools, call]
# arguments? (e.g., -o FILE)

inputs:
  outputFileType: # select type of output file: b|u|z|v
    type: string
    inputBinding:
      position: 2
      prefix: -O
      separate: false 
  ploidyType:
    type: int
    inputBinding:
      position: 3
      prefix: --ploidy
  outputFileName:
    type: string
    inputBinding:
      position: 4
      prefix: -o
  inputFileMpileup:
    type: File
    streamable: true
    inputBinding:
      position: 100
  # calling model argument
  calling_model:
    type: string
    default: "m"
    inputBinding:
      position: 5
      prefix: "-v"
      separate: false

# arguments:
  ## hard-coded parameters:
  # - prefix: -mv # output variant sites only
  #   position: 1
  # - prefix: --multiallelic-caller # alternative model for multiallelic and rare-variant calling
  #   position: 1

outputs:
  output_vcf_file:
    #type: stdout
    type: File
    outputBinding:
      glob: "$(inputs.outputFileName)"

#stdout: ERR1217015_sorted_reads.vcf
