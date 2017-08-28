cwlVersion: v1.0
class: CommandLineTool
baseCommand: node
hints:
 - class: DockerRequirement
   dockerImageId: metabric_survival
baseCommand: Rscript
arguments: [--vanilla]
inputs:
 - id: script
   type: File
   inputBinding:
    position: 1
 - id: metabric
   type: File
   inputBinding:
    position: 2
 - id: argsin
   type: File
   inputBinding:
    position: 3
outputs:
 - id: multivarTables
   type:
    type: array
    items: File
   outputBinding:
    glob: "*.tab"
 - id: km_plots
   type:
    type: array
    items: File
   outputBinding:
    glob: "*.pdf"