#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
    DockerRequirement:
        dockerPull: "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"

baseCommand: [seqtk, sample]

inputs:
    random_seed:
        type: int
        inputBinding:
            position: 1
            prefix: -s
        default: 12345
    fastq_file:
        type: File
        inputBinding:
            position: 2
    seq_size: 
        type: float
        inputBinding:
            position: 3
        default: 0.01
outputs:
    seqtk_sample_out:
        type: stdout

stdout: $(inputs.fastq_file.nameroot)