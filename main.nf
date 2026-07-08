include { VIRALCONSENSUS} from './modules/nf-core/viralconsensus/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'

process GET_REFERENCE_NAMES {
  input: path(reference)
  output: path("${reference.baseName}_names.txt")
  script:
  """
  grep ">" ${reference} | sed 's/>//g' | sed 's/ //g' > ${reference.baseName}_names.txt
  """
}

process SAMTOOLS_SUBSET {
    label 'process_low'

    conda "./modules/nf-core/samtools/view/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input: tuple val(reference_name), path(bam)
    output: tuple val("${reference_name.tokenize('|')[0]}"), val("${bam.baseName}_${reference_name.tokenize('|')[0]}"), path("${bam.baseName}_*.bam")

    script:
    def refname = reference_name.tokenize('|')[0]

    """
    samtools index ${bam}
    samtools view -h ${bam} | grep -e "${refname}" -e "^@PG" | samtools view -bS > ${bam.baseName}_${refname}.bam
    """
}

process SPLIT_REFERENCE {
    label 'process_low'
    input: tuple val(name), path(reference)
    output: tuple val(name.tokenize('|')[0]), path("${reference.baseName}_${name.tokenize('|')[0]}.fasta")

    script:
    def refname = name.tokenize('|')[0]
    """
    awk -v id="${name}" '
        /^>/ {
            keep = (\$2 == id || substr(\$1,2) == id)
        }
        keep
    ' "${reference}" > "${reference.baseName}_${refname}.fasta"
    """
}

workflow {
    main:
    // Load bam alignments
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], [file(row.bam)])
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }
    bam_ch | SAMTOOLS_INDEX
    bam_indexed_ch = bam_ch
    | join(SAMTOOLS_INDEX.out.index)

    // Load reference
    reference_ch = channel.fromPath(params.reference, checkIfExists:true)
    | map {
      n ->
      return tuple(n.baseName, n)
    }

    // Split reference by segment name
    reference_names_ch = reference_ch
    | map { n -> n.get(1) }
    | GET_REFERENCE_NAMES
    | map { n -> n.readLines()}
    | flatten

    reference_names_ch
    | combine(bam_ch | map {n -> n.get(1)})
    | SAMTOOLS_SUBSET

    reference_names_ch
    | combine(reference_ch | map { n -> n.get(1)})
    | SPLIT_REFERENCE
    | view

    input_ch = SAMTOOLS_SUBSET.out
    | join(SPLIT_REFERENCE.out)
    | map { n -> return tuple(n.get(0), n.get(1), n.get(2), n.get(3), [ ], [ ], [ ])}

    VIRALCONSENSUS(
      input_ch | map { n -> tuple(id:n.get(1), n.get(2))},
      input_ch | map { n -> tuple(n.get(0), n.get(3))},
      input_ch | map { n -> n.get(4)},
      input_ch | map { n -> n.get(5)},
      input_ch | map { n -> n.get(6)}
    )

    VIRALCONSENSUS.out.fasta | view
}
