#!/usr/bin/env nextflow


// Optional inputs implementation follows pattern
// https://github.com/nextflow-io/patterns/blob/master/optional-input.nf.
process generateStatsFigs {
    tag "${sample_id}"
    publishDir "${params.dir_out_env}/${sample_id}", \
        mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(sample_h5), path("${sample_id}.h5.*")
        each path(orf_fasta)
        each path(orf_gff)
        each path(t_rna_tsv)
        each path(codon_positions_rdata)
        each path(features_tsv)
        each path(asite_disp_length_txt)
    output:
        val sample_id, emit: finished_sample_id
        tuple val(sample_id), path("ORF_TPMs_and_counts.tsv"), emit: orf_tpms_and_counts_tsv
        tuple val(sample_id), path("metagene_start_stop_read_counts.pdf") \
            , emit: metagene_start_stop_read_counts_pdf, optional: (! params.output_pdfs)
        tuple val(sample_id), path("metagene_start_stop_read_counts.tsv") \
            , emit: metagene_start_stop_read_counts_tsv
        tuple val(sample_id), path("metagene_position_length_counts_5start.tsv") \
            , emit: metagene_position_length_counts_5start_tsv
        tuple val(sample_id), path("nt_freq_per_read_position.tsv") \
            ,emit:  nt_freq_per_read_position_tsv, optional: (! params.output_metagene_normalized_profile)
        tuple val(sample_id), path("metagene_normalized_profile_start_stop.pdf") \
            , emit: metagene_normalized_profile_start_stop_pdf, optional: (! params.output_pdfs)
        tuple val(sample_id), path("metagene_normalized_profile_start_stop.tsv") \
            , emit: metagene_normalized_profile_start_stop_tsv
        tuple val(sample_id), path("read_counts_by_length.pdf") \
            , emit: read_counts_by_length_pdf, optional: (! params.output_pdfs) 
        tuple val(sample_id), path("read_counts_by_length.tsv") \
            , emit: read_counts_by_length_tsv
        tuple val(sample_id), path("metagene_start_barplot_by_length.pdf") \
            , emit: metagene_start_barplot_by_length_pdf, optional: (! params.output_pdfs) 
        tuple val(sample_id), path("metagene_start_ribogrid_by_length.pdf") \
            , emit: metagene_start_ribogrid_by_length_pdf, optional: (! params.output_pdfs) 
        tuple val(sample_id), path("normalized_density_APEsites_per_codon.pdf") \
            ,emit: normalized_density_apesites_per_codon_pdf, optional: (! (params.is_t_rna_and_codon_positions_file && params.output_pdfs))
        tuple val(sample_id), path("normalized_density_APEsites_per_codon.tsv") \
            , emit: normalized_density_apesites_per_codon_tsv, optional: (! params.is_t_rna_and_codon_positions_file)
        tuple val(sample_id), path("normalized_density_APEsites_per_codon_long.tsv") \
            ,emit: normalized_density_apesites_per_codon_long_tsv, optional: (! params.is_t_rna_and_codon_positions_file)          
        tuple val(sample_id), path("ORF_TPMs_vs_features.pdf") \
            , emit: ORF_TPMs_vs_features_pdf, optional: (! (params.is_features_file && params.output_pdfs))
        tuple val(sample_id), path("ORF_TPMs_vs_features.tsv") \
            , emit: ORF_TPMs_vs_features_tsv, optional: (! params.is_features_file)
        tuple val(sample_id), path("read_frame_per_ORF.tsv") \
            , emit: read_frame_per_orf_tsv, optional: (! params.is_asite_disp_length_file)
        tuple val(sample_id), path("read_frame_per_ORF_filtered.tsv") \
            , emit: read_frame_per_orf_filtered_tsv, optional: (! params.is_asite_disp_length_file)
        tuple val(sample_id), path("frame_proportions_per_ORF.pdf")\
            , emit: frame_proportions_per_orf_pdf, optional: (! (params.is_asite_disp_length_file && params.output_pdfs))
    shell:
        t_rna_flag = params.is_t_rna_and_codon_positions_file \
            ? "--t-rna-file=${t_rna_tsv}" : ''
        codon_positions_flag = params.is_t_rna_and_codon_positions_file \
            ? "--codon-positions-file=${codon_positions_rdata}" : ''
        features_flag = params.is_features_file \
            ? "--features-file=${features_tsv}" : ''
        asite_disp_length_flag = params.is_asite_disp_length_file \
            ? "--asite-disp-length-file=${asite_disp_length_txt}" : ''
        count_threshold_flag = params.containsKey('count_threshold') \
            ? "--count-threshold=${params['count_threshold']}": ''
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/generate_stats_figs.R \
           --num-processes=${params.num_processes} \
           --min-read-length=${params.min_read_length} \
           --max-read-length=${params.max_read_length} \
           --buffer=${params.buffer} \
           --primary-id=${params.primary_id} \
           --dataset=${params.dataset} \
           --hd-file=${sample_h5} \
           --orf-fasta-file=${orf_fasta} \
           --output-pdfs=${params.output_pdfs} \
           --rpf=${params.rpf} \
           --output-dir=. \
           --output-metagene-normalized-profile=${params.output_metagene_normalized_profile} \
           ${t_rna_flag} \
           ${codon_positions_flag} \
           ${features_flag} \
           --orf-gff-file=${orf_gff} \
           ${asite_disp_length_flag} \
           ${count_threshold_flag}
        """
}


// Prefix sample-specific TPMs files, tpms.tsv, with sample ID so all
// sample-specific TPMs files can be staged into the same directory
// for running collateTpms.
process renameTpms {
    tag "${sample_id}"
    errorStrategy 'ignore'
    input:
        tuple val(sample_id), path(orf_tpms_and_counts_tsv)
    output:
        val(sample_id), emit: tpms_sample_id
        path("${sample_id}_tpms.tsv"), emit: tpms_sample_tsv
    shell:
        """
        cp ${orf_tpms_and_counts_tsv} ${sample_id}_tpms.tsv
        """
}

process collateTpms {
    tag "${sample_ids.join(', ')}"
    publishDir "${params.dir_out_env}", mode: 'copy', overwrite: true
    input:
        val(sample_ids)
        path(orf_tpms_and_counts_tsvs)
    output:
        path("TPMs_all_CDS_all_samples.tsv"), emit: tpms_all_cds_all_samples_tsv
        val(sample_ids), emit: collate_tpms_sample_ids
    shell:
        samples_tsvs = []
        for (i = 0; i < sample_ids.size(); i++) {
            samples_tsvs.add(sample_ids[i])
            samples_tsvs.add(orf_tpms_and_counts_tsvs[i])
        }
        samples_tsvs = samples_tsvs.join(' ')
        """
        Rscript --vanilla ${workflow.projectDir}/rscripts/collate_tpms.R \
            --tpms-file=TPMs_all_CDS_all_samples.tsv \
            ${samples_tsvs}
        """
}



process createVizParamsConfigFile {
    input:
      val(viz_params_yaml)
     output:
      path("config.yaml"), emit: viz_params_config_file_yaml
    shell:
      """
      echo "${viz_params_yaml}" > "config.yaml"
      """
}

process staticHTML {
    tag "${sample_id}"
    publishDir "${params.dir_out_env}/${sample_id}", \
    mode: 'copy', overwrite: true
    input:
      path(viz_params_config_file_yaml)
      tuple val(sample_id), \
        path(sample_metagene_start_stop_read_counts_tsv), \
        path(sample_metagene_position_length_counts_5start_tsv), \
        path(sample_read_counts_by_length_tsv), \
        path(sample_metagene_normalized_profile_start_stop_tsv), \
        path(sample_read_frame_per_orf_filtered_tsv), \
        path(sample_ORF_TPMs_vs_features_tsv), \
        path(sample_normalized_density_apesites_per_codon_long_tsv)
    output:
      val(sample_id), emit: static_html_sample_ids
      val(sample_id), emit: finished_viz_sample_id
      path("${sample_id}_output_report.html"), emit: static_html_html
  
    script:
      script = "rmarkdown::render('${workflow.projectDir}/rmarkdown/AnalysisOutputs.Rmd',"
      script += "params = list("
      script += "verbose='FALSE', "
      script += "yamlfile='\$PWD/${viz_params_config_file_yaml}', "
      script += "sampleid='!{sample_id}', "
      script += "metagene_start_stop_read_counts_data_file = '\$PWD/${sample_metagene_start_stop_read_counts_tsv}', "
      script += "metagene_position_length_counts_5start_file = '\$PWD/${sample_metagene_position_length_counts_5start_tsv}', "
      script += "read_counts_by_length_data_file='\$PWD/${sample_read_counts_by_length_tsv}', "
      script += "metagene_normalized_profile_start_stop_data_file='\$PWD/${sample_metagene_normalized_profile_start_stop_tsv}' "
      if (params.is_asite_disp_length_file) {
          script += ", read_frame_per_orf_filtered_data_file='\$PWD/${sample_read_frame_per_orf_filtered_tsv}'"
      }
      if (params.is_t_rna_and_codon_positions_file) {
          script += ", normalized_density_apesites_per_codon_long_file='\$PWD/${sample_normalized_density_apesites_per_codon_long_tsv}'"
      }
      if (params.is_features_file) {
          script += ", ORF_TPMs_vs_features_file='\$PWD/${sample_ORF_TPMs_vs_features_tsv}' "
      }
      script += "), "
      script += "intermediates_dir = '\$PWD', "
      script += "output_format = 'html_document', "
      script += "output_file = '\$PWD/${sample_id}_output_report.html')"
      """
      Rscript -e "${script}"
      """
}



// create new yaml used only for interactive visualization (riboviz/#275, riboviz/#239)
// this will write out the required params to a new yaml file
// for use by run_shiny_server.R script.
process createInteractiveVizParamsConfigFile {
    publishDir "${params.dir_out_env}", mode: 'copy', overwrite: true
    input:
      val interactive_viz_params_yaml
     output:
      path("interactive_viz_config.yaml"), emit: interactive_viz_params_config_file_yaml
    shell:
      """
      echo "${interactive_viz_params_yaml}" > "interactive_viz_config.yaml"
      """
}


process countReads {
    publishDir "${params.dir_out_env}", mode: 'copy', overwrite: true
    input:
        env PYTHONPATH
        val(ribosome_fqs_yaml)
        val(samples_ids)
    output:
        path("read_counts_per_file.tsv"), emit: read_counts_per_file_tsv
    shell:
        // 'workflow.projectDir' is directory into which outputs
        // have been published.
        // TODO: It would be preferable to:
        // 1. Stage these directories into this process's
        //    work directory. If this is possible, and how to do it,
        //    if so, has not been determined.
        // 2. Stage outputs from the processes for which reads
        //    are to be counted into this process's worif sok directory.
        //    This would require a new implementation of
        //    riboviz.tools.count_reads.
        """
        echo "${ribosome_fqs_yaml}" > ribosome_fqs.yaml
        python -m riboviz.tools.count_reads \
           -c ribosome_fqs.yaml \
           -i ${file(params.dir_in_env).toAbsolutePath()} \
           -t ${file(params.dir_tmp_env).toAbsolutePath()} \
           -o ${file(params.dir_out_env).toAbsolutePath()} \
           -r read_counts_per_file.tsv
        """
}
