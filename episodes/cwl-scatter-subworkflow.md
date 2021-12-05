## Writing Nested Workflows running in parallel

When working with workflows it is necessary to be able to apply them over multiple files (e.g., when processing a number of samples). In CWL, adding the field `ScatterFeatureRequirement` in `requirements` indicates that the workflow platform must support the `scatter` and `scatterMethod` fields that operate within the steps of the workflow. 

A "scatter" operation specifies that the associated workflow step or subworkflow should execute separately over a list of input elements. Each job making up a scatter operation is independent and may be executed concurrently. If `scatter` declares more than one input parameter, `scatterMethod` describes how to decompose the input into a discrete set of jobs. 

Notably, `scatter` is running on each step independently, but since the second step is not dependent on the first step completing all languages, we aren’t using the scatter functionality efficiently. The second step expects an array as input from the first step, so it will wait until everything in step one is finished before doing anything. 

Therefore, to use `scatter` on steps that can proceed independently of other files, we utilize the ability to nest a workflow, which includes multiple steps, as single step subworkflow within another workflow on which we are currently operating. To specify a nested workflow as part of a workflow step, `SubworkflowFeatureRequirement` must be specified in the workflow or workflow step `requirements`. It is a fatal error if a workflow directly or indirectly invokes itself as a subworkflow (recursive workflows are not allowed).

### FASTQ subsetting, compression and quality check 

In our example we are executing three command-line tools separately over two different files in FASTQ format. Our goal is to:
1. Extract a fraction (1%) of the original fastq files with seqtk.
2. Compress them (as seqtk does not automatically do so).
3. Perform quality check on them using fastqc.

***seqtk_gzip_fastqc_2.cwl***

~~~
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
~~~

### Writing the workflow

As stated above, we plan on applying multiple steps over different input files separately. So we begin by specifying `ScatterFeatureRequirement` and `SubworkflowFeatureRequirement` in the `requirements` of the workflow:

~~~
requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
~~~

This time the input will be a single array of files, defined with `type: File[]` and with the unique `id` *fastqFile_array*.

~~~
inputs:
  # seqtk input
  fastqFile_array:
    type: File[]
~~~

We start on the field `steps`, where we will specify each of the steps that will execute, by adding a unique `id` *subworkflow* for the first, and in our case the only, step. Note than in the field `run`, where the path of the `.cwl` wrapper is defined, we specify the entire nested workflow that will run as one step.

~~~
steps:
  subworkflow:
    run:
      class: Workflow
~~~

We pass one input with the `id` *fastqFile* and specify that it is a file:

~~~
inputs:
  fastqFile: 
    type: File
~~~

It is important to note how input and output `id` fields connect with each other. We pass the original workflow input *fastqFile_array* to *fastqFile* and set the *subworkflow* step to `scatter` (execute separately) over the array of files passed in *fastqFile*. Thus within the *subworkflow* step, *fastqFile* represents one of the files of *fastqFile_array*.

~~~
scatter: fastqFile
in:
  fastqFile: fastqFile_array
~~~

We start creating the steps of the nested workflow by setting the first step (*seqtk_random_selection*) and passing the `seqtk.cwl` wrapper path in `run`. The importance of carefully setting the `id` fields is also shown here. Inside the `in` field *fastq_file* is set to be *fastqFile*, which was described above. 

~~~
steps:
  # seqtk random subset selection (1%)
  seqtk_random_selection:
    run: ../cwl_tools/seqtk.cwl
    in: 
      fastq_file: fastqFile
    out: [seqtk_sample_out]
~~~

But *fastq_file* is also defined as the `id` of an input file in `seqtk.cwl` wrapper, thus connecting the input of the workflow, whose value will be set in the respective `YAML input` workflow file, with the input expected by the seqtk command-line tool in this step:

***seqtk.cwl***

~~~
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
~~~

Furthermore, note how we have set the `out` field to be *seqtk_sample_out*:

~~~
out: [seqtk_sample_out]
~~~

The selection of this `id` was not random either as it matches the output `id` of `seqtk.cwl`:

~~~
outputs:
  seqtk_sample_out:
    type: stdout
~~~

Carefully setting and using the unique values of the input and output `id` fields, allows us to connect with wrappers but also to connect different steps. This is shown here, when we created the *seqtk_gzip_step* second step and set its input to be `inputFile: seqtk_random_selection/seqtk_sample_out`. This way we directly say that the output *seqtk_sample_out* of the *seqtk_random_selection* step will be the input for *inputFile* of the *seqtk_gzip_step* step.

~~~
# gzip step for seqtk output
seqtk_gzip_step:
  run: ../cwl_tools/gzip.cwl
  in: 
    # outputFileName: seqtk_gzip_outName
    inputFile: seqtk_random_selection/seqtk_sample_out
  out: [gzip_out]
~~~

By the way, *inputFile* and *gzip_out* match the input and output `id` fields of `gzip.cwl`, respectively:

***gzip.cwl***

~~~
inputs: 
  inputFile:
    type: File
    inputBinding:
      position: 1
~~~

~~~
outputs:
  gzip_out:
    type: File
    outputBinding:
      glob: $(inputs.inputFile.basename).gz
~~~

We move on to define the final step of quality control with fastqc:

~~~
# fastqc step
  fastqc_quality_check:
    run: ../cwl_tools/bio-cwl-tools/fastqc_2.cwl
    in: 
      reads_file: seqtk_gzip_step/gzip_out
    out: [zipped_file, html_file, summary_file]
~~~

Notice here that there are more than one outputs set in the `out` field of this step. That is because fastqc produces multiple files as output and we wish to capture all of them (.zip, .html and .txt). By taking a look to the `outputs` of `fastqc_2.cwl` we see that they match those of `out` in our step and they are instructed to capture files with specific extensions based on `glob`:

***fastqc_2.cwl***

~~~
outputs:

  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
  summary_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }
~~~

Having defined the `steps` of our nested workflow, we move on with defining its `outputs`:

~~~
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
~~~

We set a unique `id` for each output, specify that we are expecting a file (`type: File`) and instruct which output (e.g., `seqtk_sample_out`) is to be expected and from which step (e.g., `seqtk_random_selection`) by using the `outputSource` field (`outputSource: seqtk_random_selection/seqtk_sample_out`).

We must not forget that we are still creating the nested workflow with the `id` *subworkflow*, which is only a step of the original workflow! As a step it contains `run`, `in`, and in this case `scatter`, which we have already defined. But we are still missing `out`. Therefore, we set `out` to be a list containing the unique `id` fields of *subworkflow* `outputs`:

~~~
out: 
  - seqtk_sample_output
  - seqtk_gzip_output
  - zipFile
  - htmlFile
  - summaryFile
~~~

Finally, the nested workflow (*subworkflow*) is complete! We can finish now our workflow by creating the final `outputs`. Once again we set unique `id` fields for each output, specify that we expect an array of files and, through the `outputSource` field, from which output of the step *subworkflow* the files will come from:

~~~
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
~~~

All that is left now is to create the input `YAML file`, containing the list of objects describing the input files, and run the CWL workflow as follows:

***seqtk_gzip_fastqc_2.yml***

~~~
fastqFile_array: [
  {class: File, path: ../raw_reads/ERR1217015_1.fastq.gz},
  {class: File, path: ../raw_reads/ERR1217015_2.fastq.gz}
]
~~~

### Executing the workflow

~~~
$ cwl-runner cwl_workflows/seqtk_gzip_fastqc_2.cwl yaml_inputs/seqtk_gzip_fastqc_2.yml
INFO /usr/bin/cwl-runner 3.0.20210124104916
INFO Resolved 'cwl_workflows/seqtk_gzip_fastqc_2.cwl' to 'file:///mnt/c/Users/CERTH/CWL/20211129/cwl_workflows/seqtk_gzip_fastqc_2.cwl'
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:185:3: Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url:
                                                                                                                       s:CreativeWork
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:193:5: Warning: checking
                                                                                                                       item
                                                                                                                       Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url:
                                                                                                                       s:PostalAddress
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:208:9: Warning: checking
                                                                                                                       item
                                                                                                                       Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url: s:Person
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:205:7: Warning: checking
                                                                                                                       item
                                                                                                                       Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url:
                                                                                                                       s:Organization
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:202:5: Warning: checking
                                                                                                                       item
                                                                                                                       Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url:
                                                                                                                       s:Organization
../../../../../../CERTH/CWL/20211129/cwl_tools/bio-cwl-tools/fastqc_2.cwl:190:3: Warning: checking
                                                                                                                       item
                                                                                                                       Warning:
                                                                                                                       Unsupported scheme
                                                                                                                       in url:
                                                                                                                       s:Organization
INFO [workflow ] start
INFO [workflow ] starting step subworkflow
INFO [step subworkflow] start
INFO [workflow subworkflow] start
INFO [workflow subworkflow] starting step seqtk_random_selection
INFO [step seqtk_random_selection] start
INFO [job seqtk_random_selection] /tmp/moec_vty$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/moec_vty,target=/sTPKdj \
    --mount=type=bind,source=/tmp/nq6uczhh,target=/tmp \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/raw_reads/ERR1217015_1.fastq.gz,target=/var/lib/cwl/stg86b35989-9016-4a31-8592-0593c6560010/ERR1217015_1.fastq.gz,readonly' \
    --workdir=/sTPKdj \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/sTPKdj \
    --cidfile=/tmp/uj69em6k/20211202161640-006664.cid \
    quay.io/biocontainers/seqtk:1.3--h5bf99c6_3 \
    seqtk \
    sample \
    -s \
    12345 \
    /var/lib/cwl/stg86b35989-9016-4a31-8592-0593c6560010/ERR1217015_1.fastq.gz \
    0.01 > /tmp/moec_vty/ERR1217015_1.fastq
INFO [job seqtk_random_selection] Max memory used: 0MiB
INFO [job seqtk_random_selection] completed success
INFO [step seqtk_random_selection] completed success
INFO [workflow subworkflow] starting step seqtk_gzip_step
INFO [step seqtk_gzip_step] start
INFO [job seqtk_gzip_step] /tmp/vyhri345$ gzip \
    -c \
    /tmp/wqgzsvut/stg2b0352a2-adee-46ee-a109-d8aa04ebe933/ERR1217015_1.fastq > /tmp/vyhri345/ERR1217015_1.fastq.gz
INFO [job seqtk_gzip_step] completed success
INFO [step seqtk_gzip_step] completed success
INFO [workflow subworkflow] starting step fastqc_quality_check
INFO [step fastqc_quality_check] start
INFO [job fastqc_quality_check] /tmp/0q7u_tco$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/0q7u_tco,target=/sTPKdj \
    --mount=type=bind,source=/tmp/twlri4tl,target=/tmp \
    --mount=type=bind,source=/tmp/vyhri345/ERR1217015_1.fastq.gz,target=/var/lib/cwl/stg3afedacb-15af-4c03-9c27-63a8bfe0bfb3/ERR1217015_1.fastq.gz,readonly \
    --workdir=/sTPKdj \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/sTPKdj \
    --cidfile=/tmp/fn0u5p6u/20211202161649-899210.cid \
    biowardrobe2/fastqc:v0.11.5 \
    fastqc \
    --extract \
    --outdir \
    . \
    /var/lib/cwl/stg3afedacb-15af-4c03-9c27-63a8bfe0bfb3/ERR1217015_1.fastq.gz
Started analysis of ERR1217015_1.fastq.gz
Approx 5% complete for ERR1217015_1.fastq.gz
Approx 15% complete for ERR1217015_1.fastq.gz
Approx 25% complete for ERR1217015_1.fastq.gz
Approx 30% complete for ERR1217015_1.fastq.gz
Approx 40% complete for ERR1217015_1.fastq.gz
Approx 50% complete for ERR1217015_1.fastq.gz
Approx 55% complete for ERR1217015_1.fastq.gz
Approx 65% complete for ERR1217015_1.fastq.gz
Approx 75% complete for ERR1217015_1.fastq.gz
Approx 80% complete for ERR1217015_1.fastq.gz
Approx 90% complete for ERR1217015_1.fastq.gz
Analysis complete for ERR1217015_1.fastq.gz
INFO [job fastqc_quality_check] Max memory used: 199MiB
INFO [job fastqc_quality_check] completed success
INFO [step fastqc_quality_check] completed success
INFO [workflow subworkflow] completed success
INFO [step subworkflow] start
INFO [workflow subworkflow_2] start
INFO [workflow subworkflow_2] starting step seqtk_random_selection_2
INFO [step seqtk_random_selection_2] start
INFO [job seqtk_random_selection_2] /tmp/i2slwxqc$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/i2slwxqc,target=/sTPKdj \
    --mount=type=bind,source=/tmp/gk3h7ics,target=/tmp \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/raw_reads/ERR1217015_2.fastq.gz,target=/var/lib/cwl/stg6b70d0b0-4a9e-403b-81b2-8c835158bd24/ERR1217015_2.fastq.gz,readonly' \
    --workdir=/sTPKdj \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/sTPKdj \
    --cidfile=/tmp/7rfqsxze/20211202161656-428241.cid \
    quay.io/biocontainers/seqtk:1.3--h5bf99c6_3 \
    seqtk \
    sample \
    -s \
    12345 \
    /var/lib/cwl/stg6b70d0b0-4a9e-403b-81b2-8c835158bd24/ERR1217015_2.fastq.gz \
    0.01 > /tmp/i2slwxqc/ERR1217015_2.fastq
INFO [job seqtk_random_selection_2] Max memory used: 0MiB
INFO [job seqtk_random_selection_2] completed success
INFO [step seqtk_random_selection_2] completed success
INFO [workflow subworkflow_2] starting step seqtk_gzip_step_2
INFO [step seqtk_gzip_step_2] start
INFO [job seqtk_gzip_step_2] /tmp/mpf8cnoz$ gzip \
    -c \
    /tmp/08oaech3/stgd43c30b5-31d4-4340-97e6-e107ef34e184/ERR1217015_2.fastq > /tmp/mpf8cnoz/ERR1217015_2.fastq.gz
INFO [job seqtk_gzip_step_2] completed success
INFO [step seqtk_gzip_step_2] completed success
INFO [workflow subworkflow_2] starting step fastqc_quality_check_2
INFO [step fastqc_quality_check_2] start
INFO [job fastqc_quality_check_2] /tmp/4lts3cey$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/4lts3cey,target=/sTPKdj \
    --mount=type=bind,source=/tmp/vjbwpavq,target=/tmp \
    --mount=type=bind,source=/tmp/mpf8cnoz/ERR1217015_2.fastq.gz,target=/var/lib/cwl/stgf50f0517-2712-4aec-8432-5128f8192cbe/ERR1217015_2.fastq.gz,readonly \
    --workdir=/sTPKdj \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/sTPKdj \
    --cidfile=/tmp/mxd9hpj8/20211202161702-633342.cid \
    biowardrobe2/fastqc:v0.11.5 \
    fastqc \
    --extract \
    --outdir \
    . \
    /var/lib/cwl/stgf50f0517-2712-4aec-8432-5128f8192cbe/ERR1217015_2.fastq.gz
Started analysis of ERR1217015_2.fastq.gz
Approx 5% complete for ERR1217015_2.fastq.gz
Approx 15% complete for ERR1217015_2.fastq.gz
Approx 25% complete for ERR1217015_2.fastq.gz
Approx 30% complete for ERR1217015_2.fastq.gz
Approx 40% complete for ERR1217015_2.fastq.gz
Approx 50% complete for ERR1217015_2.fastq.gz
Approx 55% complete for ERR1217015_2.fastq.gz
Approx 65% complete for ERR1217015_2.fastq.gz
Approx 75% complete for ERR1217015_2.fastq.gz
Approx 85% complete for ERR1217015_2.fastq.gz
Approx 90% complete for ERR1217015_2.fastq.gz
Analysis complete for ERR1217015_2.fastq.gz
INFO [job fastqc_quality_check_2] Max memory used: 171MiB
INFO [job fastqc_quality_check_2] completed success
INFO [step fastqc_quality_check_2] completed success
INFO [workflow subworkflow_2] completed success
INFO [step subworkflow] completed success
INFO [workflow ] completed success
{
    "htmlFile2": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_fastqc.html",
            "basename": "ERR1217015_1_fastqc.html",
            "class": "File",
            "checksum": "sha1$7716e4444ec780ad330db3eab9b0d8ebf140d8ad",
            "size": 289369,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_fastqc.html"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_fastqc.html",
            "basename": "ERR1217015_2_fastqc.html",
            "class": "File",
            "checksum": "sha1$7c18623105c611b20756d644e5c3d8b39854900f",
            "size": 281478,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_fastqc.html"
        }
    ],
    "seqtk_gzip_output2": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq.gz",
            "basename": "ERR1217015_1.fastq.gz",
            "class": "File",
            "checksum": "sha1$c44f71d7fa2a66beceeffccb220788242bb86694",
            "size": 908979,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq.gz"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq.gz",
            "basename": "ERR1217015_2.fastq.gz",
            "class": "File",
            "checksum": "sha1$c8ed8ec8def7a4da1b6940a9110d0573e71a772b",
            "size": 918451,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq.gz"
        }
    ],
    "seqtk_sample_output2": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq",
            "basename": "ERR1217015_1.fastq",
            "class": "File",
            "checksum": "sha1$627f39c4a3bdc2f4f11cc1b08277be0bbef501d7",
            "size": 3626489,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq",
            "basename": "ERR1217015_2.fastq",
            "class": "File",
            "checksum": "sha1$80efd304d5ce97e768e370112110f5587deb1f2f",
            "size": 3626489,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq"
        }
    ],
    "summaryFile2": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/summary.txt",
            "basename": "summary.txt",
            "class": "File",
            "checksum": "sha1$738efee3209f8bc8cf19f3024d56a47e173c8a77",
            "size": 602,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/summary.txt"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/summary.txt_2",
            "basename": "summary.txt",
            "class": "File",
            "checksum": "sha1$44943fa700dca8340f0a9adbbdfaf6415e96833c",
            "size": 602,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/summary.txt_2"
        }
    ],
    "zipFile2": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_fastqc.zip",
            "basename": "ERR1217015_1_fastqc.zip",
            "class": "File",
            "checksum": "sha1$26f4270b01aa30aa7b0c40e78388da59a98abd78",
            "size": 363023,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_fastqc.zip"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_fastqc.zip",
            "basename": "ERR1217015_2_fastqc.zip",
            "class": "File",
            "checksum": "sha1$01360fe85a6c8372d29745c162e087cafbb5c255",
            "size": 352022,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_fastqc.zip"
        }
    ]
}
INFO Final process status is success
~~~

### Resources

Detailed information are available in [CWL workflow documentation](https://www.commonwl.org/v1.2/Workflow.html). Also, CWL User Guide can be very helpful, and, especially for this example, the guides for [Nested](https://www.commonwl.org/user_guide/22-nested-workflows/index.html) and [Scattering Workflows](https://www.commonwl.org/user_guide/23-scatter-workflow/index.html).

The `fastqc_2.cwl` wrapper was adapted from [bio-cwl-tools](https://github.com/common-workflow-library/bio-cwl-tools).

FASTQ files with raw reads from DNA-Seq were retrieved from [ENA](https://www.ebi.ac.uk/ena/browser/home) (Sample Accession: SAMEA3512096).

Reference genome (FASTA) was retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/) (NCBI Reference Sequence: NC_016845.1).