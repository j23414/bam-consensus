include { VIRALCONSENSUS} from './modules/nf-core/viralconsensus/main'

process SPLIT_REFERENCE {
    tag "${meta.segment}"
    label 'process_low'
    input: tuple val(meta), path(reference)
    output: tuple val(meta), path("${reference.baseName}_*.fasta")

    script:
    def refname = meta.segment
    """
    awk -v id="${refname}" '
        /^>/ {
            keep = (\$2 ~ id || substr(\$1,2) ~ id)
        }
        keep
    ' "${reference}" > "${reference.baseName}_${refname}.fasta"
    """
}

process SAMTOOLS_SUBSET {
    tag "${meta.id},${meta.segment}"
    label 'process_low'

    conda "./modules/nf-core/samtools/view/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer']
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input: tuple val(meta), path(bam), path(reference)
    output: tuple val(meta), path("*_subset.bam"), path(reference)

    script:
    def sample = meta.id
    def refname = meta.segment

    """
    samtools index ${bam}
    samtools view -h ${bam} | grep -e "${refname}" -e "^@PG" | samtools view -bS > ${sample}_${refname}_subset.bam
    """
}

workflow {
    main:
    // Load bam alignments
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], file(row.bam))
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }

    // Load reference
    reference_ch = channel.fromPath(params.reference, checkIfExists: true)

    reference_segments_ch = reference_ch.flatMap { fasta ->
      fasta.readLines()
        .findAll { it.startsWith(">") }
        .collect { header ->
            def id = header.substring(1).split("\\|")[0]
            tuple([segment:"${id}"], fasta)
        }
    }
    | SPLIT_REFERENCE
    
    combine_ch = bam_ch
    | combine(reference_segments_ch)
    | map { sample_meta, bam, ref_meta, fasta ->
        def meta = sample_meta + ref_meta
        tuple(meta, bam, fasta)
    }
    | SAMTOOLS_SUBSET
    | map { meta, bam, fasta -> tuple([id:"${meta.id}_${meta.segment}"], bam, fasta, [], [], [])}

    VIRALCONSENSUS(
        combine_ch | map{ meta,bam,fasta,bed,save_pos,save_neg -> tuple(meta, bam)},
        combine_ch | map{ meta,bam,fasta,bed,save_pos,save_neg -> tuple(meta, file(fasta))},
        combine_ch | map{ meta,bam,fasta,bed,save_pos,save_neg -> bed},
        combine_ch | map{ meta,bam,fasta,bed,save_pos,save_neg -> save_pos},
        combine_ch | map{ meta,bam,fasta,bed,save_pos,save_neg -> save_neg}
    )
}
