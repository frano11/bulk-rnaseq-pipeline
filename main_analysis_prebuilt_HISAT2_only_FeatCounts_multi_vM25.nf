#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ----------------------
// PARAMETERS
// ----------------------
params.reads = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/nextflow_qc/results/trimmed/*.fq.gz"
params.hisat2_index = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/nextflow_alignment/genome/hisat2_index/mm10/genome"
params.annotation = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/nextflow_alignment/genome/gtf/gencode.vM25.annotation.gtf"
params.bed12 = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/nextflow_alignment/genome/BED12/gencode.vM25.annotation.bed12"
params.rrna_gtf = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/nextflow_alignment/genome/gtf/gencode.vM25_rRNA.gtf"
params.picard_path = "/Users/Frano/bin/picard.jar"
params.outdir = "results"
params.multiqc_config = "vM25_multiqc_config.yaml"

// ----------------------
// CHANNEL DEFINITIONS
// ----------------------
reads_ch = Channel
    .fromPath(params.reads)
    .map { read ->
        def id = read.getBaseName().replace('_R1_001_trimmed', '')
        tuple(id, read)
    }

index_ch      = Channel.value(params.hisat2_index)
annotation_ch = Channel.value(file(params.annotation))
bed12_ch      = Channel.value(file(params.bed12))
rrna_gtf_ch   = Channel.value(file(params.rrna_gtf))

// ----------------------
// PROCESS: HISAT2 ALIGNMENT
// ----------------------
process ALIGN_HISAT2 {

    tag "${sample_id}"
    cpus 2                    // Allocate 2 CPUs; can be increased dynamically
    maxForks 1               // Run alignments one at a time to reduce memory usage
    publishDir "${params.outdir}/Alignment", mode: 'copy', pattern: "*.sorted.bam"
    publishDir "${params.outdir}/Alignment/logs", mode: 'copy', pattern: "*.log"

    input:
    val index_prefix         // HISAT2 index path
    tuple val(sample_id), path(read)  // Sample ID and single-end FASTQ file

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),       // Sorted BAM output
          path("${sample_id}_hisat2.log")        // Alignment log

    script:
    """
    # HISAT2 alignment using the provided genome index and read file.
    # Piped into samtools to directly sort the BAM output.
    hisat2 -x ${index_prefix} \\
           -U ${read} \\
           --summary-file ${sample_id}_hisat2.log \\
           --no-softclip \\
           --pen-noncansplice 1000000 \\
           --seed 42 \\
           -p ${task.cpus} | \\
    samtools sort -@ ${task.cpus} -m 1G -o ${sample_id}.sorted.bam -
    """
}


// ----------------------
// PROCESS: REMOVE DUPLICATES (Picard MarkDuplicates)
// ----------------------
process REMOVE_DUPLICATES {
    tag "${sample_id}"
    maxForks 1
    publishDir "${params.outdir}/dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}_dedup_metrics.txt")

    script:
    """
    java -jar ${params.picard_path} MarkDuplicates \\
        I=${bam} \\
        O=${sample_id}.dedup.bam \\
        M=${sample_id}_dedup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        CREATE_INDEX=false \\
        VALIDATION_STRINGENCY=SILENT
    """
}

// ----------------------
// PROCESS: FEATURECOUNTS (gene-level quantification from deduplicated BAM)
// ----------------------
process FEATURECOUNTS {
    tag "${sample_id}"
    maxForks 1
    publishDir "${params.outdir}/FeatureCounts", mode: 'copy', pattern: "*.txt*"

    input:
    tuple val(sample_id), path(dedup_bam), path(hisat_log)
    file annotation

    output:
    tuple val(sample_id),
          path("${sample_id}_counts.txt"),
          path("${sample_id}_counts.txt.summary")

    script:
    """
    featureCounts -T 2 \\
                  -M \\
                  -s 2 \\
                  -t exon \\
                  -a ${annotation} \\
                  -o ${sample_id}_counts.txt \\
                  ${dedup_bam}
    """
}

// ----------------------
// PROCESS: QC METRICS (rRNA, strandedness, duplication)
// ----------------------
process QC_METRICS {
    tag "${sample_id}"
    maxForks 1
    publishDir "${params.outdir}/qc_metrics", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(hisat_log)
    file bed12
    file rrna_gtf

    output:
    tuple path("${sample_id}_rrna_counts.txt"),
          path("${sample_id}_strandedness.txt"),
          path("${sample_id}_raw_dup_metrics.txt")

    script:
    """
    featureCounts -T 1 -a ${rrna_gtf} -o ${sample_id}_rrna_counts.txt ${bam}

    infer_experiment.py -r ${bed12} -i ${bam} > ${sample_id}_strandedness.txt 2>&1

    java -jar ${params.picard_path} MarkDuplicates \\
        I=${bam} \\
        O=/dev/null \\
        M=${sample_id}_raw_dup_metrics.txt \\
        CREATE_INDEX=false \\
        VALIDATION_STRINGENCY=SILENT
    """
}

// ----------------------
// PROCESS: MULTIQC
// ----------------------
process MULTIQC {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path "*hisat2.log"                 // HISAT2 alignment logs
    path "*counts.txt.summary"         // featureCounts summary files
    path "*rrna_counts.txt"            // rRNA count metrics
    path "*strandedness.txt"           // Strand specificity metrics
    path "*_raw_dup_metrics.txt"       // Duplication metrics before deduplication
    path "*_dedup_metrics.txt"         // Duplication metrics after deduplication

    output:
    path "multiqc_report.html"         // Consolidated MultiQC report

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --config ${file(params.multiqc_config)} \\
        --module hisat2 \\
        --module featurecounts \\
        --module picard \\
        --module custom_content \\
        --cl-config "extra_fn_clean_exts: [ '_strandedness', '_rrna_counts', '_dedup_metrics', '_raw_dup_metrics' ]"
    """
}

// ----------------------
// WORKFLOW DEFINITION
// ----------------------
workflow {
    // Align reads using HISAT2
    aligned = ALIGN_HISAT2(index_ch, reads_ch)

    // Remove duplicates from aligned BAM files
    deduped = REMOVE_DUPLICATES(aligned.map { id, bam, log -> tuple(id, bam) })

    // Count features using featureCounts
    counted = FEATURECOUNTS(
        deduped.combine(aligned.map { id, bam, log -> tuple(id, log) }),
        annotation_ch
    )

    // Generate QC metrics
    qc_metrics = QC_METRICS(aligned, bed12_ch, rrna_gtf_ch)

    // Extract specific outputs for MultiQC
    hisat2_logs     = aligned.map { it[2] }        // HISAT2 logs
    counts_summary  = counted.map { it[2] }        // featureCounts summaries
    rrna_counts     = qc_metrics.map { it[0] }     // rRNA counts
    strand_metrics  = qc_metrics.map { it[1] }     // Strand specificity metrics
    raw_dup_metrics = qc_metrics.map { it[2] }     // Pre-deduplication metrics
    dedup_metrics   = deduped.map { it[2] }        // Post-deduplication metrics

    // Run MultiQC with all collected metrics
    MULTIQC(
        hisat2_logs,
        counts_summary,
        rrna_counts,
        strand_metrics,
        raw_dup_metrics,
        dedup_metrics
    )
}

