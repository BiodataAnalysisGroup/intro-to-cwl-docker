## Preparing secondary files before running bioinformatic analyses.

### Constructing FM-index for the reference genome

Often before running a workflow in bioinformatics it is necessary to prepare various files (e.g., indexes) which will be used during the analysis by one or more of the envoked command-line tools of the workflow. 

In this example we are running `bwa index` through CWL to prepare the required secondary files before mapping raw reads to the reference genome with `bwa mem`.

First we prepare the CWL wrapper for the command-line tool bwa index, called `BWA-Index.cwl`, or modify an existing wrapper. Commonly used tools have been implemented in CWL and are readily available from various sources, like the [bio-cwl-tools](https://github.com/common-workflow-library/bio-cwl-tools) repository, from where we adopted many of the tools used in these examples.

***BWA-Index.cwl***

~~~
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
~~~

We specify the version of CWL and that we are running a command-line tool:

~~~
cwlVersion: v1.0
class: CommandLineTool
~~~

For simplification and compatibility we specify that we are using Docker containers rather than pre-installed tools. This is done by using `DockerRequirement` in the `requirements` or `hints` fields, and in this case the container `quay.io/biocontainers/bwa:0.7.17--ha92aebf_3` through `dockerPull`.

~~~
requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
~~~

We define the command that we want using the `baseCommand` field as well as the `inputs` that we are passing to our tool. The two inputs used in this example require defining a unique input `id` for each, which in this case are *InputFile* and *IndexName*, for reference genome and the name of the output file, respectively. 

~~~
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
~~~

Note that *InputFile* and *IndexName* are defined as `type: File` and `type: string`, as we are passing a FASTA file and a string to them as values, respectively. 

Also, we can, optionally, define the expected format(s) of the input file(s) in the `format` field, which in this case is FASTA. If we do, specifying `$namespaces` and `$schemas` is required:

~~~
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
~~~

With the field `inputBinding` we can add prefixes, like `-p`, and define the order of the arguments, like the `position: 200` which is a simple way to set the input reference genome to be preceded by any other input with lower `position`.

The input *algoType* is a good example of an optional input, for which a value is not necessary to be defined latter in the input `YAML` file. CWL recognizes this if we: 
    - add `- "null"` to the `type` field as the example, or
    - add a `?` after the value of `type` (e.g., `type: File?` or `type: string?`)

~~~
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
~~~

Finally, we must define the expected `outputs`:

~~~
outputs:
  index:
    doc: all secondary files required for bwa mem
    type:
      type: array
      items: File
    outputBinding:
      glob: "$(inputs.IndexName)*"
~~~

As with `inputs`, a unique `id` must be specified for each output, which in our case is *index*. Note that we are expecting multiple files as output and, therefore, we have to:
    - specify the field `type` to be an array of files
    - use `glob` from the `outputBinding` files to capture multiple files that have the string from the *IndexName* input in their file names.
    
For the latter, we can access directly the *IndexName* string value with `$(inputs.IndexName)`, and specify that we are expecting different extensions by adding `*`.

As our CWL is ready before we can execute it, we have to create an input `YAML` file, where we tell the tool what values to expect for its `inputs`:

***BWA-Index.yml***

~~~
# bwa index input
InputFile: 
  class: File
  format: edam:format_1929  # FASTA
  path: ../NC_016845.fasta
IndexName: NC_016845.fasta
~~~

We use the unique input `id` to pass the respective values for each input. For `InputFile`, we specify again that this is a file, with `class`, and that it is of FASTA format, necessary if the respective input has defined the `format` field. Note here that while `format` can be a list that defines multiple compatible file formats in the `.cwl` wrapper (e.g., FASTQ or FASTA), it can only take one value, that is the format of the input file (e.g., FASTA) in the input `YAML` file.

Additionaly, the `path` for the FASTA file is set for `InputFile`, while `IndexName` only takes a string value as it was defined to be of `type: string` in *BWA-Index.cwl*.

By invoking `cwl-runner` we can run the tool like this and get the expected output files:

~~~
$cwl-runner cwl_tools/bio-cwl-tools/BWA-Index.cwl yaml_inputs/BWA-Index.yml
INFO /usr/bin/cwl-runner 3.0.20210124104916
INFO Resolved 'cwl_tools/bio-cwl-tools/BWA-Index.cwl' to 'file:///mnt/c/Users/CWL/20211129/cwl_tools/bio-cwl-tools/BWA-Index.cwl'
INFO [job BWA-Index.cwl] /tmp/wp3q6n0i$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/wp3q6n0i,target=/rgDuaw \
    --mount=type=bind,source=/tmp/hyw5oq2f,target=/tmp \
    '--mount=type=bind,source=/mnt/c/Users/CWL/20211129/NC_016845.fasta,target=/var/lib/cwl/stg90d8a3de-f290-4cae-b86b-9a0c19748f1a/NC_016845.fasta,readonly' \    --workdir=/rgDuaw \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/rgDuaw \
    --cidfile=/tmp/482n3owq/20211201123553-872048.cid \
    quay.io/biocontainers/bwa:0.7.17--ha92aebf_3 \
    bwa \
    index \
    -p \
    NC_016845.fasta \
    /var/lib/cwl/stg90d8a3de-f290-4cae-b86b-9a0c19748f1a/NC_016845.fasta
[bwa_index] Pack FASTA... 0.09 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 4.23 seconds elapse.
[bwa_index] Update BWT... 0.07 sec
[bwa_index] Pack forward-only FASTA... 0.05 sec
[bwa_index] Construct SA from BWT and Occ... 1.03 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -p NC_016845.fasta /var/lib/cwl/stg90d8a3de-f290-4cae-b86b-9a0c19748f1a/NC_016845.fasta
[main] Real time: 5.829 sec; CPU: 5.556 sec
INFO [job BWA-Index.cwl] Max memory used: 0MiB
INFO [job BWA-Index.cwl] completed success
{
    "index": [
        {
            "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.amb",
            "basename": "NC_016845.fasta.amb",
            "class": "File",
            "checksum": "sha1$20c0e67325d4e0f9ba0483eb9356317e6b4e8d10",
            "size": 24,
            "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.amb"
        },
        {
            "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.ann",
            "basename": "NC_016845.fasta.ann",
            "class": "File",
            "checksum": "sha1$3eb9f3ebce8a13341b99304dec580d03a52a26eb",
            "size": 115,
            "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.ann"
        },
        {
            "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.bwt",
            "basename": "NC_016845.fasta.bwt",
            "class": "File",
            "checksum": "sha1$cec4937a5ebe2d2c119cc0da2041b2d770cade3c",
            "size": 5334020,
            "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.bwt"
        },
        {
            "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.pac",
            "basename": "NC_016845.fasta.pac",
            "class": "File",
            "checksum": "sha1$81bd78382d78b76a435003bc6e4d7a9c9e2cd458",
            "size": 1333487,
            "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.pac"
        },
        {
            "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.sa",
            "basename": "NC_016845.fasta.sa",
            "class": "File",
            "checksum": "sha1$dc8f9c9957af590186eea741401746127a382b79",
            "size": 2667024,
            "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.sa"
        }
    ]
}
INFO Final process status is success
$
$ ls -1 NC_016845.fasta.*
NC_016845.fasta.amb
NC_016845.fasta.ann
NC_016845.fasta.bwt
NC_016845.fasta.fai
NC_016845.fasta.pac
NC_016845.fasta.sa 
~~~

### Index reference sequence in FASTA format

In the example workflow (`cwl-variant-calling-example-workflow.md`), an additional index `.fai` of the reference sequence is required before running our analysis, and it can be produced with the `samtools faidx` tool.

As with `BWA-Index.cwl` we create or modify an existing wrapper (in this case adopted from bio-cwl-tools):

***samtools_faidx.cwl***

~~~
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]

baseCommand: [ samtools, faidx ]

inputs:
  sequences:
    type: File
    doc: Input FASTA file
    format: edam:format_1929

arguments:
   - $(inputs.sequences.basename)

outputs:
  sequences_with_index:
    type: File
    format: $(inputs.sequences.format)
    secondaryFiles:
     - .fai
    outputBinding:
      glob: $(inputs.sequences.basename)
  sequences_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename).fai

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
~~~

Note the usage of `InitialWorkDirRequirement` in `requirements`, where files and/or subfolders, necessary to set before running the tool, can be listed.

~~~
requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]
~~~

Interestingly, the required FASTA file is defined in `inputs` with a unique `id`, and the latter is used by several fields within the tool (`$(inputs.sequences)`) to access its attributes (`basename` and `format`). However, `inputBinding` is not defined but rather the file is passed within the command by using its attribute `$(inputs.sequences.basename)` in the `arguments` field.

~~~
arguments:
   - $(inputs.sequences.basename)
~~~

Furthermore, although optional, the input FASTA file is defined as the expected output *sequences_with_index* with its necessary `secondaryFile` being the second expected output *sequences_index*:

~~~
outputs:
  sequences_with_index:
    type: File
    format: $(inputs.sequences.format)
    secondaryFiles:
     - .fai
    outputBinding:
      glob: $(inputs.sequences.basename)
  sequences_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename).fai
~~~

Similarly, a `YAML_input` file is created and the tool is executed with `cwl-runner`:

***samtools_faidx.yml***

~~~
sequences: 
  class: File
  format: edam:format_1929  # FASTA
  path: ../NC_016845.fasta
~~~

~~~
$cwl-runner cwl_tools/bio-cwl-tools/samtools_faidx.cwl yaml_inputs/samtools_faidx.yml
INFO /usr/bin/cwl-runner 3.0.20210124104916
INFO Resolved 'cwl_tools/bio-cwl-tools/samtools_faidx.cwl' to 'file:///mnt/c/Users/CWL/20211129/cwl_tools/bio-cwl-tools/samtools_faidx.cwl'
INFO [job samtools_faidx.cwl] /tmp/vky2awf3$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/vky2awf3,target=/BZKWqG \
    --mount=type=bind,source=/tmp/qbwmx07o,target=/tmp \
    '--mount=type=bind,source=/mnt/c/Users/CWL/20211129/NC_016845.fasta,target=/BZKWqG/NC_016845.fasta,readonly' \
    --workdir=/BZKWqG \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/BZKWqG \
    --cidfile=/tmp/fpk441mm/20211201133347-174331.cid \
    kerstenbreuer/samtools:1.7 \
    samtools \
    faidx \
    NC_016845.fasta
INFO [job samtools_faidx.cwl] Max memory used: 0MiB
INFO [job samtools_faidx.cwl] completed success
{
    "sequences_index": {
        "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.fai",
        "basename": "NC_016845.fasta.fai",
        "class": "File",
        "checksum": "sha1$fb42b7b6000c7beceac88aa701cd053b7e1af9e3",
        "size": 29,
        "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.fai"
    },
    "sequences_with_index": {
        "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta",
        "basename": "NC_016845.fasta",
        "class": "File",
        "checksum": "sha1$a844314d0bde6e1e18d19052a0b7a6a5f1e7c4ed",
        "size": 5410232,
        "secondaryFiles": [
            {
                "basename": "NC_016845.fasta.fai",
                "location": "file:///mnt/c/Users/CWL/20211129/NC_016845.fasta.fai",
                "class": "File",
                "checksum": "sha1$fb42b7b6000c7beceac88aa701cd053b7e1af9e3",
                "size": 29,
                "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta.fai"
            }
        ],
        "format": "http://edamontology.org/format_1929",
        "path": "/mnt/c/Users/CWL/20211129/NC_016845.fasta"
    }
}
INFO Final process status is success
$
$ ls -1 | grep ".fai"
NC_016845.fasta.fai
~~~
