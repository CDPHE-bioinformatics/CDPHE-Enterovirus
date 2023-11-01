version 1.0

workflow EVD68_seq_results {

    input {
        Array[String] sample_name
        Array[File?] cov_out
        Array[File?] percent_cvg_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[File] terra_data_table_path_array
        
        # python scripts
        File evd68_concat_seq_results_py

    }

    # secret variables - for static values convert from array to single entity
    String project_name = select_all(project_name_array)[0]
    File terra_data_table_path = select_all(terra_data_table_path_array)[0]
    String out_dir = select_all(out_dir_array)[0]

    call results_table {
      input:
        sample_name = sample_name,
        evd68_concat_seq_results_py = evd68_concat_seq_results_py,
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        terra_data_table_path = terra_data_table_path
    }

    call transfer {
      input:
          out_dir = out_dir,
          sequencing_results_csv = results_table.sequencing_results_csv
    }

    output {
        File sequencing_results_csv = results_table.sequencing_results_csv

    }
}

task results_table {

    input {
      Array[String] sample_name
      File evd68_concat_seq_results_py
      Array[File] cov_out
      Array[File] percent_cvg_csv
      String project_name
      File terra_data_table_path

    }

    command <<<
    python ~{evd68_concat_seq_results_py} \
        --sample_name_array "~{write_lines(sample_name)}" \
        --cov_out_files "~{write_lines(cov_out)}" \
        --percent_cvg_files "~{write_lines(percent_cvg_csv)}" \
        --project_name "~{project_name}" \
        --terra_data_table_path "~{terra_data_table_path}" \

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"

    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
}

task transfer {
    input {
        String out_dir
        File sequencing_results_csv
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/

    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}

