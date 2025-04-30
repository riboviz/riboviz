#!/usr/bin/env nextflow

/*
Sample file (fq_files)-specific processes.
*/

process cutAdapters {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_fq) 
    output:
        tuple val(sample_id), file("trim.fq"), emit: cut_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o trim.fq ${sample_fq} -j 0
        """
}

process extractUmis {
    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "${dir_tmp}/${sample_id}", \
        mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(sample_id), file(sample_fq) \
            from cut_fq_branch.umi_fq
    output:
        tuple val(sample_id), file("extract_trim.fq"), emit: umi_extract_fq
    shell:
        """
        umi_tools extract -I ${sample_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S extract_trim.fq
        """
}

/*
Multiplexed files (multiplex_fq_files)-specific processes.
*/

process cutAdaptersMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${dir_tmp}", mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_fq) \
            from multiplex_id_fq.collect{ id, file -> [id, file] }
    output:
        tuple val(multiplex_id), file("${multiplex_id}_trim.fq"), emit: cut_multiplex_fq
    shell:
        """
        cutadapt --trim-n -O 1 -m 5 -a ${params.adapters} \
            -o ${multiplex_id}_trim.fq ${multiplex_fq} -j 0
        """
}


process extractUmisMultiplex {
    tag "${multiplex_id}"
    errorStrategy 'ignore'
    publishDir "${dir_tmp}", mode: publish_index_tmp_type, overwrite: true
    input:
        tuple val(multiplex_id), file(multiplex_fq) \
            from cut_multiplex_fq_branch.umi_fq
    output:
        tuple val(multiplex_id), file("${multiplex_id}_extract_trim.fq"), emit: umi_extract_multiplex_fq
    when:
        params.extract_umis && is_multiplexed
    shell:
        """
        umi_tools extract -I ${multiplex_fq} \
            --bc-pattern="${params.umi_regexp}" \
            --extract-method=regex -S ${multiplex_id}_extract_trim.fq
        """
}


process demultiplex {
    tag "${multiplex_id}"
    publishDir "${dir_tmp}/${multiplex_id}_deplex", \
        mode: publish_index_tmp_type, overwrite: true
    errorStrategy 'ignore'
    input:
        // Use '.toString' to prevent changing hashes of
        // 'workflow.projectDir' triggering reexecution of this
        // process if 'nextflow run' is run with '-resume'.
        env PYTHONPATH from workflow.projectDir.toString()
        tuple val(multiplex_id), file(multiplex_fq)
        each file(sample_sheet_tsv)
    output:
        tuple val(multiplex_id), file("num_reads.tsv"), emit: demultiplex_num_reads_tsv
        file("*.f*"), emit: demultiplex_fq
    shell:
        """
        python -m riboviz.tools.demultiplex_fastq \
            -1 ${multiplex_fq} -s ${sample_sheet_tsv} -o . -m 2
        """
}