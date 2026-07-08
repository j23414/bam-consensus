# bam-consensus

A modular workflow for generating the consensus sequence from alignment information.

## Usage

```
nextflow run j23414/bam-consensus \
  --bam [path/*.bam] \
  --samplesheet [path/bam_samplesheet.csv] \
  --reference [path/reference.fasta] \
  --outdir "consensus-results" \
  -profile stjude
```
