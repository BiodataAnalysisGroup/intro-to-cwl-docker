#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

requirements:
  ScatterFeatureRequirement: {}

inputs:
  message_array: string[] 

steps:
  echo:
    run: echo.cwl
    scatter: message
    in:
      message: message_array
    out: []

outputs: []