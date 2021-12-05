## Basic variant calling workflow

In this example we create a basic variant calling CWL workflow, which serves as an example on how to connect different tools to analyze diffent biological data in different formats and meaningful information.

### Trimming raw reads, mapping and variant calling

In our example we are executing several command-line tools, taking as initial input a random fraction (1%) of two paired-end FASTQ files (DNA-Seq), in order to:
1. Perform trimming and quality control of the initial raw reads.
2. Map the trimmed reads to a reference genome.
3. Convert SAM to BAM format.
4. Sort BAM file by coordinates.
5. Perform variant calling with bcftools.

***variant_calling_example_workflow.cwl***

~~~
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
# input of trim_galore step:  
  raw_fastq_file_1: 
    type: File
  raw_fastq_file_2: 
    type: File
# input for bwa mem
  ref_genome_2nd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
# input for samtools view sam2bam
  sort_type:
    type: boolean
    default: false
# input for bcftools mpileup
  ref_genome_3rd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .fai
  output_file_type:
    type: string
    default: "u"
# input for bcftools call
  output_file_type_2nd:
    type: string
    default: "v"
  type_of_ploidy: 
    type: int
  output_file_name:
    type: string
  type_of_model:
    type: string
    default: "m"

outputs: 
# output of trim_galore step:
  fastq_1_trimmed_file:
    type: File
    outputSource: trimming_reads/fastq1_trimmed
  fastq_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed
  fastq_unpaired_1_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq1_trimmed_unpaired
  fastq_unpaired_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed_unpaired
  trim_galoreLog:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trim_galore_log
  trimmed_fastqc_html_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_html
  trimmed_fastqc_zip_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_zip
# output for bwa mem
  bwa_mem_file:
    type: File
    outputSource: bwa_mem_mapping_reads/reads_stdout
# output for samtools view sam2bam
  sam_to_bam_file:
    type: File
    outputSource: sam_to_bam_conversion/bam
# output for samtools sort
  bam_sorted_file:
    type: File
    outputSource: bam_file_sorting/bam_sorted
# output for bcftools mpileup call pipeline
  bcftools_vcf_file:
    type: File
    outputSource: bcftools_call_step/output_vcf_file


steps:
# trim_galore trimming of reads
  trimming_reads:
    run: ../cwl_tools/bio-cwl-tools/trim_galore.cwl
    in:
      fastq1: raw_fastq_file_1
      fastq2: raw_fastq_file_2
    out:
      - fastq1_trimmed 
      - fastq2_trimmed
      - fastq1_trimmed_unpaired 
      - fastq2_trimmed_unpaired 
      - trim_galore_log
      - trimmed_fastqc_html
      - trimmed_fastqc_zip
# bwa mem mapping reads:
  bwa_mem_mapping_reads:
    run: ../cwl_tools/bio-cwl-tools/BWA-Mem.cwl
    in: 
      InputFile_1: trimming_reads/fastq1_trimmed # trimmed_reads_1
      InputFile_2: trimming_reads/fastq2_trimmed
      Index: ref_genome_2nd
    out: [reads_stdout]
# SAM_to_BAM_conversion:
  sam_to_bam_conversion:
    run: ../cwl_tools/bio-cwl-tools/samtools_view_sam2bam.cwl
    in:
      sam: bwa_mem_mapping_reads/reads_stdout # sam_file
    out: [bam]
# BAM_sorting:
  bam_file_sorting:
    run: ../cwl_tools/bio-cwl-tools/samtools_sort.cwl
    in:
      bam_unsorted: sam_to_bam_conversion/bam  #bam_file
      by_name: sort_type
    out: [bam_sorted]
# bcftools_mpileup_step:
  bcftools_mpileup_step:
    run: ../cwl_tools/bcftools_mpileup.cwl
    in:
      refGenome: ref_genome_3rd
      outputFileType: output_file_type
      inputBam: bam_file_sorting/bam_sorted # bam_file_2nd
    out: [pileupOut]
# bcftools_call_step:
  bcftools_call_step:
    run: ../cwl_tools/bcftools_call.cwl
    in:
      outputFileType: output_file_type_2nd
      ploidyType: type_of_ploidy
      outputFileName: output_file_name
      inputFileMpileup: bcftools_mpileup_step/pileupOut
      calling_model: type_of_model
    out: [output_vcf_file]

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
~~~

### Writing the workflow

We begin by setting CWL version and specifying that we are writing a workflow:

~~~
cwlVersion: v1.0
class: Workflow
~~~

The field `inputs` is subsequently filled with the input(s) that the tool of each step will require. The specific values of these inputs will be latter set in the input `YAML file`.

Note the functionality offered by certain fields, including:
- `secondaryFile`, that provides a pattern or expression specifying files or directories that should be included alongside the primary file. Secondary files may be required or optional.
- `format`, which must be one or more IRIs of concept nodes that represents file formats which are allowed as input to this parameter, preferrably defined within an ontology. If no ontology is available, file formats may be tested by exact match.
- `default`, which specifies the default value to use for this parameter if the parameter is missing from the input object, or if the value of the parameter in the input object is *null*. Default values are applied before evaluating expressions.

~~~
inputs:
# input of trim_galore step:  
  raw_fastq_file_1: 
    type: File
  raw_fastq_file_2: 
    type: File
# input for bwa mem
  ref_genome_2nd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
# input for samtools view sam2bam
  sort_type:
    type: boolean
    default: false
# input for bcftools mpileup
  ref_genome_3rd: # secondary files required
    type: File
    format: edam:format_1929 # FASTA
    secondaryFiles: 
      - .fai
  output_file_type:
    type: string
    default: "u"
# input for bcftools call
  output_file_type_2nd:
    type: string
    default: "v"
  type_of_ploidy: 
    type: int
  output_file_name:
    type: string
  type_of_model:
    type: string
    default: "m"
~~~

We then proceed to write the `steps` of this workflow. At the first step, we set to run the wrapper `trim_galore.cwl` to trim and perform quality check on the paired-end fastq files. The path of the tool is specified in `run`, inputs in `in` and the list of outputs in `out`:

~~~
# trim_galore trimming of reads
trimming_reads:
run: ../cwl_tools/bio-cwl-tools/trim_galore.cwl
in:
    fastq1: raw_fastq_file_1
    fastq2: raw_fastq_file_2
out:
    - fastq1_trimmed 
    - fastq2_trimmed
    - fastq1_trimmed_unpaired 
    - fastq2_trimmed_unpaired 
    - trim_galore_log
    - trimmed_fastqc_html
    - trimmed_fastqc_zip
~~~

Note that the `id` unique input:

~~~
in:
    fastq1: raw_fastq_file_1
    fastq2: raw_fastq_file_2
~~~

and output values: 

~~~
out:
    - fastq1_trimmed 
    - fastq2_trimmed
    - fastq1_trimmed_unpaired 
    - fastq2_trimmed_unpaired 
    - trim_galore_log
    - trimmed_fastqc_html
    - trimmed_fastqc_zip
~~~

specified here are not set in random but match those of the invoked `trim_galore.cwl` wrapper:

***trim_galore.cwl***

`inputs`:

~~~
inputs:
  # main input
  fastq1:
    doc: |
      raw reads in fastq format; can be gzipped;
      if paired end, the file contains the first reads;
      if single end, the file contains all reads
    type: File
    inputBinding:
      position: 10
  fastq2:
    doc: |
      (optional) raw reads in fastq format; can be gzipped;
      if paired end, the file contains the second reads;
      if single end, the file does not exist
    type: File?
    inputBinding:
      position: 11
~~~

`outputs`:

~~~
outputs:
  fastq1_trimmed:
    type: File
    outputBinding:
      glob: |
        ${
            if ( inputs.fastq2 == null  ){ return "*trimmed.fq*" }
            else { return "*val_1.fq*" }
        }
  fastq2_trimmed:    
    type: File?
    outputBinding:
      glob: "*val_2.fq*"
  fastq1_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: "*unpaired_1.fq*"
  fastq2_trimmed_unpaired:    
    type: File?
    outputBinding:
      glob: "*unpaired_2.fq*"
  trim_galore_log: # can be used by multiqc
    type:
      type: array # since one or two matches (single/paired end)
      items: File
    outputBinding:
      glob:  "*trimming_report.txt"
  trimmed_fastqc_html:
    doc: html report of post-trimming fastqc
    type:
      type: array # since one or two matches (single/paired end)
      items: File
    outputBinding:
      glob: "*fastqc.html"
  trimmed_fastqc_zip:
    doc: all data of post-trimming fastqc e.g. figures
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputBinding:
      glob: "*fastqc.zip"
~~~

Carefully setting and using the unique values of the input and output `id` fields, allows us to connect with wrappers but also to connect different steps. The latter can be seen in the next step where we set two of the three inputs (*InputFile_1* and *InputFile_2*) to take the values of two of the direct outputs of the previous step:

~~~
# bwa mem mapping reads:
bwa_mem_mapping_reads:
run: ../cwl_tools/bio-cwl-tools/BWA-Mem.cwl
in: 
    InputFile_1: trimming_reads/fastq1_trimmed # trimmed_reads_1
    InputFile_2: trimming_reads/fastq2_trimmed
    Index: ref_genome_2nd
out: [reads_stdout]
~~~

Similarly, we specify the next two steps of converting SAM to BAM format and sorting the BAM file by coordinates:

~~~
# SAM_to_BAM_conversion:
  sam_to_bam_conversion:
    run: ../cwl_tools/bio-cwl-tools/samtools_view_sam2bam.cwl
    in:
      sam: bwa_mem_mapping_reads/reads_stdout # sam_file
    out: [bam]
# BAM_sorting:
  bam_file_sorting:
    run: ../cwl_tools/bio-cwl-tools/samtools_sort.cwl
    in:
      bam_unsorted: sam_to_bam_conversion/bam  #bam_file
      by_name: sort_type
    out: [bam_sorted]
~~~

Notably, the output of the `bcftools_mpileup.cwl` wrapper, which is called in the step *bcftools_mpileup_step*, is passed to the standard output stream and captured by specifying the field `stdout`:

***bcftools_mpileup.cwl***

`outputs`:

~~~
outputs:
  pileupOut:
    type: stdout
~~~

We utilize the output `id` value *pileupOut* to capture and pass this output directly as input to the next step *bcftools_call_step*:

~~~
# bcftools_mpileup_step:
  bcftools_mpileup_step:
    run: ../cwl_tools/bcftools_mpileup.cwl
    in:
      refGenome: ref_genome_3rd
      outputFileType: output_file_type
      inputBam: bam_file_sorting/bam_sorted # bam_file_2nd
    out: [pileupOut]
# bcftools_call_step:
  bcftools_call_step:
    run: ../cwl_tools/bcftools_call.cwl
    in:
      outputFileType: output_file_type_2nd
      ploidyType: type_of_ploidy
      outputFileName: output_file_name
      inputFileMpileup: bcftools_mpileup_step/pileupOut
      calling_model: type_of_model
    out: [output_vcf_file]
~~~

Producing this way the variant calling pipeline of *bcftools mpileup | bcftools call*.

Having concluded with writing the `steps` of this workflow we can proceed with capturing the `outputs` of each step:

~~~
outputs: 
# output of trim_galore step:
  fastq_1_trimmed_file:
    type: File
    outputSource: trimming_reads/fastq1_trimmed
  fastq_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed
  fastq_unpaired_1_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq1_trimmed_unpaired
  fastq_unpaired_2_trimmed_file:
    type: File?
    outputSource: trimming_reads/fastq2_trimmed_unpaired
  trim_galoreLog:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trim_galore_log
  trimmed_fastqc_html_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_html
  trimmed_fastqc_zip_report:
    type:
      type: array
      items: File # since one or two matches (single/paired end)
    outputSource: trimming_reads/trimmed_fastqc_zip
# output for bwa mem
  bwa_mem_file:
    type: File
    outputSource: bwa_mem_mapping_reads/reads_stdout
# output for samtools view sam2bam
  sam_to_bam_file:
    type: File
    outputSource: sam_to_bam_conversion/bam
# output for samtools sort
  bam_sorted_file:
    type: File
    outputSource: bam_file_sorting/bam_sorted
# output for bcftools mpileup call pipeline
  bcftools_vcf_file:
    type: File
    outputSource: bcftools_call_step/output_vcf_file
~~~

Interestingly, we define the `type` of the expected output to be:
- a file, `type: File`
- or an array of items that are files:
~~~
type:
  type: array
  items: File
~~~

With tools like *Trim_galore* where one or more files are expected as input (single- or paired-end file(s)), and thus the respective number of files is produced as output, we can set the an input or an output are optional by adding `?` directly after `type` (e.g., `type: File?`).

Note that again the unique `id` values come into play to connect this time different `steps` and `outputs`. For example the first output that we are trying to capture is the first fastq file that has been trimmed by *Trim_galore*. Therefore, we give a unique `id` to the first output (*fastq_1_trimmed_file*), set that we are expecting a file (`type: File`), and use the field `outputSource` to specify that from the step *trimming_reads* we are capturing the output *fastq1_trimmed* (`outputSource: trimming_reads/fastq1_trimmed`):

~~~
# output of trim_galore step:
  fastq_1_trimmed_file:
    type: File
    outputSource: trimming_reads/fastq1_trimmed
~~~

Now that are workflow is complete, we only need to set the input values in the input `YAML file`:

~~~
# trim_galore input
raw_fastq_file_1:
  class: File
  path: ../ERR1217015_1.fastq.gz
raw_fastq_file_2:
  class: File
  path: ../ERR1217015_2.fastq.gz
# bwa mem input
ref_genome_2nd:
  class: File
  path: ../NC_016845.fasta
  format: edam:format_1929 # FASTA
# bcftools mpileup input
ref_genome_3rd:
    class: File
    path: ../NC_016845.fasta
    format: edam:format_1929 # FASTA
# output_file_type: "u"
# bcftools call input
output_file_type_2nd: "v"
type_of_ploidy: 1
#type_of_model: "m"
output_file_name: "ERR1217015_sorted_reads.vcf"
~~~

and run the workflow as follows:

### Executing the workflow

~~~
cwl-runner cwl_workflows/variant_calling_example_workflow.cwl yaml_inputs/variant_calling_example_workflow.yml
INFO /usr/bin/cwl-runner 3.0.20210124104916
INFO Resolved 'cwl_workflows/variant_calling_example_workflow.cwl' to 'file:///mnt/c/Users/CERTH/CWL/20211129/cwl_workflows/variant_calling_example_workflow.cwl'  
Could not load extension schema http://edamontology.org/EDAM_1.18.owl: Error fetching http://edamontology.org/EDAM_1.18.owl: HTTPConnectionPool(host='edamontology.org', port=80): Max retries exceeded with url: /EDAM_1.18.owl (Caused by NewConnectionError('<urllib3.connection.HTTPConnection object at 0x7ff47395f3d0>: Failed to establish a new connection: [Errno -3] Temporary failure in name resolution'))
Warning: Field `$schemas` contains undefined reference to `http://edamontology.org/EDAM_1.18.owl`
Could not load extension schema http://edamontology.org/EDAM_1.18.owl: Error fetching http://edamontology.org/EDAM_1.18.owl: HTTPConnectionPool(host='edamontology.org', port=80): Max retries exceeded with url: /EDAM_1.18.owl (Caused by NewConnectionError('<urllib3.connection.HTTPConnection object at 0x7ff4736a3910>: Failed to establish a new connection: [Errno -3] Temporary failure in name resolution'))
Warning: Field `$schemas` contains undefined reference to `http://edamontology.org/EDAM_1.18.owl`
WARNING Workflow checker warning:
../../../../../../CERTH/CWL/20211129/cwl_workflows/variant_calling_example_workflow.cwl:102:9: Source
                                                                                                                                     'fastq2_trimmed' of
                                                                                                                                     type ["null",      
                                                                                                                                     "File"] may be     
                                                                                                                                     incompatible       
../../../../../../CERTH/CWL/20211129/cwl_workflows/variant_calling_example_workflow.cwl:113:7:   with sink        
                                                                                                                                       'InputFile_2' of 
                                                                                                                                       type "File"      
INFO [workflow ] start
INFO [workflow ] starting step trimming_reads
INFO [step trimming_reads] start
WARNING [job trimming_reads] Skipping Docker software container '--memory' limit despite presence of ResourceRequirement with ramMin and/or ramMax setting. Consider running with --strict-memory-limit for increased portability assurance.
INFO [job trimming_reads] /tmp/v35xijrt$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/v35xijrt,target=/wFvDJH \
    --mount=type=bind,source=/tmp/f3dpbccu,target=/tmp \   
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq.gz,target=/var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz,readonly' \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq.gz,target=/var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz,readonly' \
    --workdir=/wFvDJH \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/gpd7eu3n/20211205125522-436420.cid \
    kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7 \
    trim_galore \
    --fastqc_args \
    '"--noextract"' \
    --paired \
    --retain_unpaired \
    --length_1 \
    35 \
    --length_2 \
    35 \
    --stringency \
    1 \
    --length \
    20 \
    --quality \
    20 \
    /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz \
    /var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz
No quality encoding type selected. Assuming that the data provided uses Sanger encoded Phred scores (default)

Path to Cutadapt set as: 'cutadapt' (default)
1.14
Cutadapt seems to be working fine (tested command 'cutadapt --version')


AUTO-DETECTING ADAPTER TYPE
===========================
Attempting to auto-detect adapter type from the first 1 million sequences of the first file (>> /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz <<)

Found perfect matches for the following adapter sequences:
Adapter type    Count   Sequence        Sequences analysed      Percentage
Illumina        0       AGATCGGAAGAGC   11791   0.00
Nextera 0       CTGTCTCTTATA    11791   0.00
smallRNA        0       TGGAATTCTCGG    11791   0.00
Unable to auto-detect most prominent adapter from the first specified file (count Illumina: 0, count Nextera: 0)
Defaulting to Illumina universal adapter ( AGATCGGAAGAGC ). Specify -a SEQUENCE to avoid this behavior).

Writing report to 'ERR1217015_1.fastq.gz_trimming_report.txt'

SUMMARISING RUN PARAMETERS
==========================
Input filename: /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz
Trimming mode: paired-end
Trim Galore version: 0.4.4
Cutadapt version: 1.14
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; default (inconclusive auto-detection))
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
Length cut-off for read 1: 35 bp (default)
Length cut-off for read 2: 35 bb (default)
Running FastQC on the data once trimming has completed
Running FastQC with the following extra arguments: '"--noextract"'
Output file(s) will be GZIP compressed

Writing final adapter and quality trimmed output to ERR1217015_1_trimmed.fq.gz


  >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz <<<
This is cutadapt 1.14 with Python 2.7.12
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 0.36 s (31 us/read; 1.96 M reads/minute).

=== Summary ===

Total reads processed:                  11,791
Reads with adapters:                     2,979 (25.3%)
Reads written (passing filters):        11,791 (100.0%)

Total basepairs processed:     1,473,875 bp
Quality-trimmed:                   5,285 bp (0.4%)
Total written (filtered):      1,464,297 bp (99.4%)

=== Adapter 1 ===

Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 2979 times.

No. of allowed errors:
0-9 bp: 0; 10-13 bp: 1

Bases preceding removed adapters:
  A: 28.0%
  C: 34.0%
  G: 22.9%
  T: 15.1%
  none/other: 0.0%

Overview of removed sequences
length  count   expect  max.err error counts
1       2170    2947.8  0       2170
2       641     736.9   0       641
3       121     184.2   0       121
4       30      46.1    0       30
5       9       11.5    0       9
6       5       2.9     0       5
84      2       0.0     1       0 2
115     1       0.0     1       0 1


RUN STATISTICS FOR INPUT FILE: /var/lib/cwl/stgae38a0d0-a4c7-44da-a7b3-28d244146e71/ERR1217015_1.fastq.gz
=============================================
11791 sequences processed in total
The length threshold of paired-end sequences gets evaluated later on (in the validation step)

Writing report to 'ERR1217015_2.fastq.gz_trimming_report.txt'

SUMMARISING RUN PARAMETERS
==========================
Input filename: /var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz
Trimming mode: paired-end
Trim Galore version: 0.4.4
Cutadapt version: 1.14
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC' (Illumina TruSeq, Sanger iPCR; default (inconclusive auto-detection))
Maximum trimming error rate: 0.1 (default)
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
Length cut-off for read 1: 35 bp (default)
Length cut-off for read 2: 35 bb (default)
Running FastQC on the data once trimming has completed
Running FastQC with the following extra arguments: '"--noextract"'
Output file(s) will be GZIP compressed

Writing final adapter and quality trimmed output to ERR1217015_2_trimmed.fq.gz


  >>> Now performing quality (cutoff 20) and adapter trimming in a single pass for the adapter sequence: 'AGATCGGAAGAGC' from file /var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz <<<
This is cutadapt 1.14 with Python 2.7.12
Command line parameters: -f fastq -e 0.1 -q 20 -O 1 -a AGATCGGAAGAGC /var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz
Trimming 1 adapter with at most 10.0% errors in single-end mode ...
Finished in 0.37 s (32 us/read; 1.90 M reads/minute).

=== Summary ===

Total reads processed:                  11,791
Reads with adapters:                     2,843 (24.1%)
Reads written (passing filters):        11,791 (100.0%)

Total basepairs processed:     1,473,875 bp
Quality-trimmed:                  10,662 bp (0.7%)
Total written (filtered):      1,459,244 bp (99.0%)

=== Adapter 1 ===

Sequence: AGATCGGAAGAGC; Type: regular 3'; Length: 13; Trimmed: 2843 times.

No. of allowed errors:
0-9 bp: 0; 10-13 bp: 1

Bases preceding removed adapters:
  A: 30.5%
  C: 32.3%
  G: 21.2%
  T: 16.0%
  none/other: 0.0%

Overview of removed sequences
length  count   expect  max.err error counts
1       2034    2947.8  0       2034
2       640     736.9   0       640
3       119     184.2   0       119
4       34      46.1    0       34
5       13      11.5    0       13
6       1       2.9     0       1
9       1       0.0     0       1
82      1       0.0     1       0 1


RUN STATISTICS FOR INPUT FILE: /var/lib/cwl/stg4f5dce4a-e4cb-4576-8ac9-d67c16305a53/ERR1217015_2.fastq.gz
=============================================
11791 sequences processed in total
The length threshold of paired-end sequences gets evaluated later on (in the validation step)

Validate paired-end files ERR1217015_1_trimmed.fq.gz and ERR1217015_2_trimmed.fq.gz
file_1: ERR1217015_1_trimmed.fq.gz, file_2: ERR1217015_2_trimmed.fq.gz


>>>>> Now validing the length of the 2 paired-end infiles: ERR1217015_1_trimmed.fq.gz and ERR1217015_2_trimmed.fq.gz <<<<<
Writing validated paired-end read 1 reads to ERR1217015_1_val_1.fq.gz
Writing validated paired-end read 2 reads to ERR1217015_2_val_2.fq.gz

Writing unpaired read 1 reads to ERR1217015_1_unpaired_1.fq.gz
Writing unpaired read 2 reads to ERR1217015_2_unpaired_2.fq.gz

Total number of sequences analysed: 11791

Number of sequence pairs removed because at least one read was shorter than the length cutoff (20 bp): 17 (0.14%)


  >>> Now running FastQC on the validated data ERR1217015_1_val_1.fq.gz<<<

Started analysis of ERR1217015_1_val_1.fq.gz
Approx 5% complete for ERR1217015_1_val_1.fq.gz
Approx 15% complete for ERR1217015_1_val_1.fq.gz
Approx 25% complete for ERR1217015_1_val_1.fq.gz
Approx 30% complete for ERR1217015_1_val_1.fq.gz
Approx 40% complete for ERR1217015_1_val_1.fq.gz
Approx 50% complete for ERR1217015_1_val_1.fq.gz
Approx 55% complete for ERR1217015_1_val_1.fq.gz
Approx 65% complete for ERR1217015_1_val_1.fq.gz
Approx 75% complete for ERR1217015_1_val_1.fq.gz
Approx 80% complete for ERR1217015_1_val_1.fq.gz
Approx 90% complete for ERR1217015_1_val_1.fq.gz
Analysis complete for ERR1217015_1_val_1.fq.gz

  >>> Now running FastQC on the validated data ERR1217015_2_val_2.fq.gz<<<

Started analysis of ERR1217015_2_val_2.fq.gz
Approx 5% complete for ERR1217015_2_val_2.fq.gz
Approx 15% complete for ERR1217015_2_val_2.fq.gz
Approx 25% complete for ERR1217015_2_val_2.fq.gz
Approx 30% complete for ERR1217015_2_val_2.fq.gz
Approx 40% complete for ERR1217015_2_val_2.fq.gz
Approx 50% complete for ERR1217015_2_val_2.fq.gz
Approx 55% complete for ERR1217015_2_val_2.fq.gz
Approx 65% complete for ERR1217015_2_val_2.fq.gz
Approx 75% complete for ERR1217015_2_val_2.fq.gz
Approx 85% complete for ERR1217015_2_val_2.fq.gz
Approx 90% complete for ERR1217015_2_val_2.fq.gz
Analysis complete for ERR1217015_2_val_2.fq.gz
Deleting both intermediate output files ERR1217015_1_trimmed.fq.gz and ERR1217015_2_trimmed.fq.gz

====================================================================================================

INFO [job trimming_reads] Max memory used: 211MiB
INFO [job trimming_reads] completed success
INFO [step trimming_reads] completed success
INFO [workflow ] starting step bwa_mem_mapping_reads
INFO [step bwa_mem_mapping_reads] start
INFO [job bwa_mem_mapping_reads] /tmp/5cyuic83$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/5cyuic83,target=/wFvDJH \
    --mount=type=bind,source=/tmp/josnn7d6,target=/tmp \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta,readonly' \    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.amb,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta.amb,readonly' \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.ann,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta.ann,readonly' \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.bwt,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta.bwt,readonly' \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.pac,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta.pac,readonly' \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.sa,target=/var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta.sa,readonly' \
    --mount=type=bind,source=/tmp/v35xijrt/ERR1217015_1_val_1.fq.gz,target=/var/lib/cwl/stga540f031-2ddd-40c6-8fd2-4721704aa822/ERR1217015_1_val_1.fq.gz,readonly \
    --mount=type=bind,source=/tmp/v35xijrt/ERR1217015_2_val_2.fq.gz,target=/var/lib/cwl/stge9f40081-b62d-4437-a90e-eaa990e0b0b7/ERR1217015_2_val_2.fq.gz,readonly \
    --workdir=/wFvDJH \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/5ma6igwt/20211205125555-315469.cid \
    quay.io/biocontainers/bwa:0.7.17--ha92aebf_3 \
    bwa \
    mem \
    /var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta \
    /var/lib/cwl/stga540f031-2ddd-40c6-8fd2-4721704aa822/ERR1217015_1_val_1.fq.gz \
    /var/lib/cwl/stge9f40081-b62d-4437-a90e-eaa990e0b0b7/ERR1217015_2_val_2.fq.gz > /tmp/5cyuic83/unsorted_reads.sam
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 23548 sequences (2921510 bp)...
[M::mem_pestat] # candidate unique pairs for (FF, FR, RF, RR): (0, 10609, 2, 0)
[M::mem_pestat] skip orientation FF as there are not enough pairs
[M::mem_pestat] analyzing insert size distribution for orientation FR...
[M::mem_pestat] (25, 50, 75) percentile: (446, 498, 561)
[M::mem_pestat] low and high boundaries for computing mean and std.dev: (216, 791)
[M::mem_pestat] mean and std.dev: (506.27, 84.81)
[M::mem_pestat] low and high boundaries for proper pairs: (101, 906)
[M::mem_pestat] skip orientation RF as there are not enough pairs
[M::mem_pestat] skip orientation RR as there are not enough pairs
[M::mem_process_seqs] Processed 23548 reads in 1.793 CPU sec, 1.837 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem /var/lib/cwl/stg8236d7db-536c-4e21-8ac6-39af52dd84e8/NC_016845.fasta /var/lib/cwl/stga540f031-2ddd-40c6-8fd2-4721704aa822/ERR1217015_1_val_1.fq.gz /var/lib/cwl/stge9f40081-b62d-4437-a90e-eaa990e0b0b7/ERR1217015_2_val_2.fq.gz
[main] Real time: 2.197 sec; CPU: 2.041 sec
INFO [job bwa_mem_mapping_reads] Max memory used: 25MiB
INFO [job bwa_mem_mapping_reads] completed success
INFO [step bwa_mem_mapping_reads] completed success
INFO [workflow ] starting step sam_to_bam_conversion
INFO [step sam_to_bam_conversion] start
WARNING [job sam_to_bam_conversion] Skipping Docker software container '--memory' limit despite presence of ResourceRequirement with ramMin and/or ramMax setting. Consider running with --strict-memory-limit for increased portability assurance.
INFO [job sam_to_bam_conversion] /tmp/r4fjvdzy$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/r4fjvdzy,target=/wFvDJH \
    --mount=type=bind,source=/tmp/pxu10_ck,target=/tmp \
    --mount=type=bind,source=/tmp/5cyuic83/unsorted_reads.sam,target=/var/lib/cwl/stg5e881461-75fe-4bb6-85be-eafcc35d9448/unsorted_reads.sam,readonly \
    --workdir=/wFvDJH \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/6tr9nhc1/20211205125559-840899.cid \
    kerstenbreuer/samtools:1.7 \
    samtools \
    view \
    -h \
    -b \
    /var/lib/cwl/stg5e881461-75fe-4bb6-85be-eafcc35d9448/unsorted_reads.sam > /tmp/r4fjvdzy/unsorted_reads.bam
INFO [job sam_to_bam_conversion] Max memory used: 1MiB
INFO [job sam_to_bam_conversion] completed success
INFO [step sam_to_bam_conversion] completed success
INFO [workflow ] starting step bam_file_sorting
INFO [step bam_file_sorting] start
WARNING [job bam_file_sorting] Skipping Docker software container '--memory' limit despite presence of ResourceRequirement with ramMin and/or ramMax setting. Consider running with --strict-memory-limit for increased portability assurance.
INFO [job bam_file_sorting] /tmp/_ca4jhtb$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/_ca4jhtb,target=/wFvDJH \
    --mount=type=bind,source=/tmp/7t924g_8,target=/tmp \
    --mount=type=bind,source=/tmp/r4fjvdzy/unsorted_reads.bam,target=/var/lib/cwl/stg7a91ee8b-6ab7-4bf7-b9ef-4081da70de4d/unsorted_reads.bam,readonly \
    --workdir=/wFvDJH \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/va23hk9y/20211205125602-238744.cid \
    kerstenbreuer/samtools:1.7 \
    samtools \
    sort \
    /var/lib/cwl/stg7a91ee8b-6ab7-4bf7-b9ef-4081da70de4d/unsorted_reads.bam > /tmp/_ca4jhtb/sorted_reads.bam
INFO [job bam_file_sorting] Max memory used: 10MiB
INFO [job bam_file_sorting] completed success
INFO [step bam_file_sorting] completed success
INFO [workflow ] starting step bcftools_mpileup_step
INFO [step bcftools_mpileup_step] start
INFO [job bcftools_mpileup_step] /tmp/lp61mz_o$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/lp61mz_o,target=/wFvDJH \
    --mount=type=bind,source=/tmp/rgqbr9g4,target=/tmp \
    --mount=type=bind,source=/tmp/_ca4jhtb/sorted_reads.bam,target=/var/lib/cwl/stgfd950eec-ba2a-491e-9695-d92ec6d82d6f/sorted_reads.bam,readonly \
    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta,target=/var/lib/cwl/stgc3427f76-0b3f-446e-94d8-a197a87b41a8/NC_016845.fasta,readonly' \    '--mount=type=bind,source=/mnt/c/Users/CERTH/CWL/20211129/NC_016845.fasta.fai,target=/var/lib/cwl/stgc3427f76-0b3f-446e-94d8-a197a87b41a8/NC_016845.fasta.fai,readonly' \
    --workdir=/wFvDJH \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/s6qniwr9/20211205125604-661660.cid \
    quay.io/biocontainers/bcftools:1.14--h88f3f91_0 \
    bcftools \
    mpileup \
    -f \
    /var/lib/cwl/stgc3427f76-0b3f-446e-94d8-a197a87b41a8/NC_016845.fasta \
    -Ou \
    /var/lib/cwl/stgfd950eec-ba2a-491e-9695-d92ec6d82d6f/sorted_reads.bam > /tmp/lp61mz_o/906c58bc207039a4848033fb8a150aedfed86537
[mpileup] 1 samples in 1 input files
[mpileup] maximum number of reads per input file set to -d 250
INFO [job bcftools_mpileup_step] Max memory used: 42MiB
INFO [job bcftools_mpileup_step] completed success
INFO [step bcftools_mpileup_step] completed success
INFO [workflow ] starting step bcftools_call_step
INFO [step bcftools_call_step] start
INFO [job bcftools_call_step] /tmp/a1j4a47v$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/a1j4a47v,target=/wFvDJH \
    --mount=type=bind,source=/tmp/srjae5hd,target=/tmp \
    --mount=type=bind,source=/tmp/lp61mz_o/906c58bc207039a4848033fb8a150aedfed86537,target=/var/lib/cwl/stge90b52a2-2ea3-48d8-95cf-3c9ba32c2b54/906c58bc207039a4848033fb8a150aedfed86537,readonly \      
    --workdir=/wFvDJH \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/wFvDJH \
    --cidfile=/tmp/zyx0a0tx/20211205125618-578415.cid \
    quay.io/biocontainers/bcftools:1.14--h88f3f91_0 \
    bcftools \
    call \
    -Ov \
    --ploidy \
    1 \
    -o \
    ERR1217015_sorted_reads.vcf \
    -vm \
    /var/lib/cwl/stge90b52a2-2ea3-48d8-95cf-3c9ba32c2b54/906c58bc207039a4848033fb8a150aedfed86537
INFO [job bcftools_call_step] Max memory used: 0MiB
INFO [job bcftools_call_step] completed success
INFO [step bcftools_call_step] completed success
INFO [workflow ] completed success
{
    "bam_sorted_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/sorted_reads.bam",
        "basename": "sorted_reads.bam",
        "class": "File",
        "checksum": "sha1$c8e7c555c7caf4a45c83ff23672ec720ba87254d",
        "size": 1921444,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/sorted_reads.bam"
    },
    "bcftools_vcf_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_sorted_reads.vcf",
        "basename": "ERR1217015_sorted_reads.vcf",
        "class": "File",
        "checksum": "sha1$03db32b9cd869e1dc2510f7924cb37f21d5d83e5",
        "size": 192304,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_sorted_reads.vcf"
    },
    "bwa_mem_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/unsorted_reads.sam",
        "basename": "unsorted_reads.sam",
        "class": "File",
        "checksum": "sha1$743057db491e87ced5b81e644859d6bd46bc4679",
        "size": 8320970,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/unsorted_reads.sam"
    },
    "fastq_1_trimmed_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1.fq.gz",
        "basename": "ERR1217015_1_val_1.fq.gz",
        "class": "File",
        "checksum": "sha1$207e7996913b8399254399a82cb7fa692ae38cf3",
        "size": 900461,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1.fq.gz"
    },
    "fastq_2_trimmed_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2.fq.gz",
        "basename": "ERR1217015_2_val_2.fq.gz",
        "class": "File",
        "checksum": "sha1$14c7be5ea678c766d1b3ce416867735bd23622ea",
        "size": 907981,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2.fq.gz"
    },
    "fastq_unpaired_1_trimmed_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_unpaired_1.fq.gz",
        "basename": "ERR1217015_1_unpaired_1.fq.gz",
        "class": "File",
        "checksum": "sha1$3b6d6999d6cd8329601e2524d6684451c26c8d5f",
        "size": 1562,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_unpaired_1.fq.gz"
    },
    "fastq_unpaired_2_trimmed_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_unpaired_2.fq.gz",
        "basename": "ERR1217015_2_unpaired_2.fq.gz",
        "class": "File",
        "checksum": "sha1$7faac864cce90269cf69e9b8e048666bcc096bfa",
        "size": 171,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_unpaired_2.fq.gz"
    },
    "sam_to_bam_file": {
        "location": "file:///mnt/c/Users/CERTH/CWL/20211129/unsorted_reads.bam",
        "basename": "unsorted_reads.bam",
        "class": "File",
        "checksum": "sha1$4dec2838e1ee46634bfe4c98a7bc199faf84a02f",
        "size": 1932536,
        "path": "/mnt/c/Users/CERTH/CWL/20211129/unsorted_reads.bam"
    },
    "trim_galoreLog": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq.gz_trimming_report.txt",
            "basename": "ERR1217015_1.fastq.gz_trimming_report.txt",
            "class": "File",
            "checksum": "sha1$536784ec126c1caf56ac2243f4098c57d8c8ae9b",
            "size": 2090,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1.fastq.gz_trimming_report.txt"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq.gz_trimming_report.txt",
            "basename": "ERR1217015_2.fastq.gz_trimming_report.txt",
            "class": "File",
            "checksum": "sha1$87f3733a42c33d7ffe6e21fae15598a974b6f242",
            "size": 2285,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2.fastq.gz_trimming_report.txt"
        }
    ],
    "trimmed_fastqc_html_report": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1_fastqc.html",
            "basename": "ERR1217015_1_val_1_fastqc.html",
            "class": "File",
            "checksum": "sha1$89c95646096939aaefde304e1c0ce6b21542812a",
            "size": 249352,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1_fastqc.html"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2_fastqc.html",
            "basename": "ERR1217015_2_val_2_fastqc.html",
            "class": "File",
            "checksum": "sha1$c651dd4ad65577e8400eeb05bd307d28208a2f25",
            "size": 250323,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2_fastqc.html"
        }
    ],
    "trimmed_fastqc_zip_report": [
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1_fastqc.zip",
            "basename": "ERR1217015_1_val_1_fastqc.zip",
            "class": "File",
            "checksum": "sha1$afa8a05db3cdd111a8ec4a39feeb6fa792f27dbc",
            "size": 316347,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_1_val_1_fastqc.zip"
        },
        {
            "location": "file:///mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2_fastqc.zip",
            "basename": "ERR1217015_2_val_2_fastqc.zip",
            "class": "File",
            "checksum": "sha1$bff85716d0e736e552ec161f35db9657fdab19d2",
            "size": 317392,
            "path": "/mnt/c/Users/CERTH/CWL/20211129/ERR1217015_2_val_2_fastqc.zip"
        }
    ]
}
INFO Final process status is success
~~~

### Resources

Detailed information are available in [CWL workflow documentation](https://www.commonwl.org/v1.2/Workflow.html). Also, CWL User Guide can be very helpful, and, especially for this example, the guide for [Writing Workfloes](https://www.commonwl.org/user_guide/21-1st-workflow/index.html).

Several of the wrappers used in this example were adapted from [bio-cwl-tools](https://github.com/common-workflow-library/bio-cwl-tools).

FASTQ files with raw reads from DNA-Seq were retrieved from [ENA](https://www.ebi.ac.uk/ena/browser/home) (Sample Accession: SAMEA3512096).

Reference genome (FASTA) was retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/) (NCBI Reference Sequence: NC_016845.1).