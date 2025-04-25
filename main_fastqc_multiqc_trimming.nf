#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ========================
// Parameters & Directories
// ========================
params.input_dir = "/Users/Frano/Desktop/Bioinfo_2025/250127_Doppelganger/April_2025_Bulk_RNAseq/fastq_gz"

params.results_dir          = "./results"
params.fastqc_raw_dir       = "${params.results_dir}/fastqc_raw"
params.trimmed_dir          = "${params.results_dir}/trimmed"
params.trimmed_reports_dir  = "${params.trimmed_dir}/reports"
params.fastqc_trimmed_dir   = "${params.trimmed_dir}/fastqc_trimmed"
params.multiqc_raw_dir      = "${params.results_dir}/multiqc_raw"
params.multiqc_trimmed_dir  = "${params.results_dir}/multiqc_trimmed"
params.nextflow_reports_dir = "${params.results_dir}/nextflow_reports"

// ========================
// Input Channels
// ========================
Channel.fromPath("${params.input_dir}/*.fastq.gz")
    .map { file ->
        def md5File = file.parent.resolve("${file.name}.md5")
        tuple(file.simpleName, file, md5File)
    }
    .set { all_fastq }

Channel.fromPath("${params.input_dir}/*.fastq.gz.md5")
    .map { md5File ->
        def fastqFile = md5File.parent.resolve(md5File.name - ".md5")
        if (!fastqFile.exists()) tuple(md5File.simpleName, md5File)
    }
    .filter { it != null }
    .set { orphaned_md5 }

// ========================
// Processes
// ========================
process check_md5 {
    tag { base }
    errorStrategy 'ignore'
    cpus 1
    memory '1.GB'
    time '10.min'

    input:
    tuple val(base), path(fastq), path(md5)

    output:
    path fastq, emit: validated
    path "${base}_md5.log", emit: log

    script:
    """
    set +e
    log_file="${base}_md5.log"

    if [ -f "${md5}" ]; then
        if command -v md5sum >/dev/null 2>&1; then
            computed=\$(md5sum "${fastq}" | cut -d ' ' -f1)
        elif command -v md5 >/dev/null 2>&1; then
            computed=\$(md5 -q "${fastq}")
        else
            echo "ERROR: No MD5 tool found" > "\$log_file"
            cp -f "${fastq}" .
            exit 0
        fi

        expected=\$(awk '{print \$1}' "${md5}")

        if [ "\$expected" = "\$computed" ]; then
            echo "SUCCESS: ${base} verified" > "\$log_file"
        else
            echo "ERROR: MD5 mismatch for ${base}" > "\$log_file"
        fi
    else
        echo "WARNING: Missing MD5 for ${base}" > "\$log_file"
    fi

    cp -f "${fastq}" .
    exit 0
    """
}

process check_orphaned_md5 {
    tag { base }
    errorStrategy 'ignore'
    cpus 1
    memory '1.GB'
    time '5.min'

    input:
    tuple val(base), path(md5)

    output:
    path "${base}_orphaned.log", emit: log

    script:
    """
    echo "ERROR: Orphaned MD5 for ${base}" > "${base}_orphaned.log"
    exit 0
    """
}

process fastqc_raw {
    cpus 2
    memory '2.GB'
    time '30.min'
    publishDir "${params.fastqc_raw_dir}", mode: 'copy'

    input:
    path fastq

    output:
    path "${fastq.simpleName}_fastqc.zip", emit: reports
    path "${fastq.simpleName}_fastqc.html"

    script:
    """
    fastqc -q "${fastq}"
    """
}

process trim_galore {
    cpus 4
    memory '8.GB'
    time '1.h'
    tag { fastq.simpleName }

    publishDir "${params.trimmed_dir}", pattern: "*_trimmed.fq.gz", mode: 'copy'
    publishDir "${params.trimmed_reports_dir}", pattern: "*trimming_report.txt", mode: 'copy'

    input:
    path fastq

    output:
    tuple val(fastq.simpleName), path("*_trimmed.fq.gz"), path("*trimming_report.txt"), emit: trimmed

    script:
    """
    trim_galore \\
      --quality 20 \\
      --length 25 \\
      --stringency 5 \\
      --illumina \\
      --clip_R1 10 \\
      --three_prime_clip_R1 1 \\
      --cores ${task.cpus} \\
      --output_dir . \\
      "${fastq}"
    """
}

process fastqc_trimmed {
    cpus 2
    memory '2.GB'
    time '30.min'
    publishDir "${params.fastqc_trimmed_dir}", mode: 'copy'

    input:
    tuple val(base), path(trimmed_fastq), path(trimming_report)

    output:
    path "${trimmed_fastq.simpleName}_fastqc.zip", emit: reports
    path "${trimmed_fastq.simpleName}_fastqc.html"

    script:
    """
    fastqc -q "${trimmed_fastq}"
    """
}

process multiqc_raw {
    cpus 2
    memory '4.GB'
    time '20.min'
    publishDir "${params.multiqc_raw_dir}", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report_raw.html"

    script:
    """
    mkdir -p raw_reports
    cp -r ${reports.join(' ')} raw_reports/
    multiqc -f raw_reports -n multiqc_report_raw.html -o .
    """
}

process multiqc_trimmed {
    cpus 2
    memory '4.GB'
    time '20.min'
    publishDir "${params.multiqc_trimmed_dir}", mode: 'copy'

    input:
    path reports

    output:
    path "multiqc_report_trimmed.html"

    script:
    """
    mkdir -p trimmed_reports
    cp -r ${reports.join(' ')} trimmed_reports/
    multiqc -f trimmed_reports -n multiqc_report_trimmed.html -o .
    """
}

process combine_md5_logs {
    publishDir "${params.results_dir}", mode: 'copy'

    input:
    path logs

    output:
    path "verify_md5_integrity.txt"

    script:
    """
    echo "MD5 Verification Report" > verify_md5_integrity.txt
    echo "Generated: \$(date '+%Y-%m-%d %H:%M:%S')" >> verify_md5_integrity.txt
    echo "----------------------------------------" >> verify_md5_integrity.txt
    cat *.log >> verify_md5_integrity.txt 2>/dev/null || echo "No validation logs found." >> verify_md5_integrity.txt
    """
}

// ========================
// Workflow
// ========================
workflow {
    check_md5(all_fastq)
    check_orphaned_md5(orphaned_md5)

    fastqc_raw(check_md5.out.validated)

    trim_galore(check_md5.out.validated)
    fastqc_trimmed(trim_galore.out.trimmed)

    multiqc_raw(
        fastqc_raw.out.reports.collect()
    )

    multiqc_trimmed(
        fastqc_trimmed.out.reports
            .mix(trim_galore.out.trimmed.collect { it[2] })
            .collect()
    )

    combine_md5_logs(
        check_md5.out.log
            .mix(check_orphaned_md5.out.log)
            .collect()
    )
}
