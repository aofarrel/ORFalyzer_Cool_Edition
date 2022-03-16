class: CommandLineTool
id: "ORFalyzer"
label: "ORFalyzer"
cwlVersion: v1.1
doc: |
    Identify ORFs and the hypothetical properties of amino acids.

dct:creator:
  foaf:name: "Ash O'Farrell"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/aofarrel/orfalyzer:python3.9.10slim"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4092
    outdirMin: 512000

inputs:
  mem_gb:
    type: File
    default: 4
    doc: "Input file"
    inputBinding:
      position: 1

  bam_input:
    type: string
    doc: "Name of the output file"
    inputBinding:
      position: 2

outputs:
  output_file:
    type: File
    outputBinding:
      glob: bamstats_report.zip
    doc: "A text file reporting ORFs and their properties as amino acids."


baseCommand: ["python3", "/orfalyzer/orfalyzer.py"]