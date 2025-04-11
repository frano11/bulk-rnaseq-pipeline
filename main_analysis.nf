#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ========================
// Configuration Parameters
// ========================
params {
    genome_fasta   = "${projectDir}/genome/mouse.fa"
    gtf            = "${projectDir}/genome/mouse.gtf"
    trimmed_dir    = "${projectDir}/results/trimmed"
    genome_index   = "${projectDir}/genome/star_index"
    alignment_dir  = "${projectDir}/results/alignment"
}

// ========================
// Channels
// ========================
// Channel for trimmed FASTQ files
Channel.fromPath("${params.trimmed_dir}/*_trimmed.fq.gz")
    .map { file -> 
        def sample = file.simpleName - ~/_trimmed$/
        tuple(sample, file)
    }
    .set { trimmed_samples }

// ========================
// Processes
// ========================
process star_genome_index {
    tag "STAR_Genome_Index"
    label 'high_mem'
    publishDir params.genome_index, mode: 'copy'

    input:
    path genome_fasta
    path gtf

    output:
    path "${params.genome_index}", emit: index

    script:
    """
    mkdir -p ${params.genome_index}
    STAR --runMode genomeGenerate \\
         --genomeDir ${params.genome_index} \\
         --genomeFastaFiles ${genome_fasta} \\
         --sjdbGTFfile ${gtf} \\
         --sjdbOverhang 74 \\
         --runThreadN ${task.cpus}
    """
}

process star_align {
    tag { sample }
    publishDir params.alignment_dir, mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".Log.final.out")) "logs/$filename"
        else filename
    }

    input:
    tuple val(sample), path(trimmed_fq)
    path genome_index

    output:
    path "${sample}.*", emit: alignment

    script:
    """
    STAR --genomeDir ${genome_index} \\
         --readFilesIn ${trimmed_fq} \\
         --readFilesCommand zcat \\
         --outSAMtype BAM SortedByCoordinate \\
         --quantMode GeneCounts \\
         --outFileNamePrefix ${sample}. \\
         --runThreadN ${task.cpus}
    """
}

// ========================
// Workflow
// ========================
workflow {
    // Build genome index
    star_index = star_genome_index(params.genome_fasta, params.gtf)

    // Align trimmed reads
    star_align(trimmed_samples, star_index.out.index)
}