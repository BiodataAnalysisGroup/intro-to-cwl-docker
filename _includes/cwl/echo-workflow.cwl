#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

inputs:
  message: string 

steps:
  echo:
    run: echo.cwl
    in:
      message: message
    out: []

outputs: []