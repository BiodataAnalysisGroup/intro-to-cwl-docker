#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gzip, -c]

requirements:
  InlineJavascriptRequirement: {}

inputs: 
  inputFile:
    type: File
    inputBinding:
      position: 1
  # outputFileName: # to be able to run for pipe input...
  #   type: string
    
stdout: $(inputs.inputFile.basename).gz # $(inputs.outputFileName).gz

outputs:
  gzip_out:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.basename).gz # $(inputs.outputFileName).gz
