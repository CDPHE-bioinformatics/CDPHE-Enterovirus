version 1.0

workflow EVD68_seq_results {

    input {
        Array[String] sample_name
        Array[File?] cov_out
        Array[File?] percent_cvg_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        # Array[String?] assembler_version_array
        # Array[File] workbook_path_array
        Array[File] terra_data_table_path_array
        Array[File] assembly_software_file_array
        
        # python scripts
        File ev-d68_seq_results_py

    }

    # secret variables - for static values convert from array to single entity
    String project_name = select_all(project_name_array)[0]
    File terra_data_table_path = select_all(terra_data_table_path_array)[0]
    String out_dir = select_all(out_dir_array)[0]
    File assembly_software_file = select_all(assembly_software_file_array)[0]
    # String assembler_version = select_all(assembler_version_array)[0]
    # File workbook_path = select_all(workbook_path_array)[0]

    call results_table {
      input:
        sample_name = sample_name,
        ev-d68_seq_results_py = ev-d68_seq_results_py,
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        terra_data_table_path = terra_data_table_path,
        assembly_software_file = assembly_software_file
    }

    call transfer {
      input:
          out_dir = out_dir,
          sequencing_results_csv = results_table.sequencing_results_csv,
          wgs_horizon_report_csv = results_table.wgs_horizon_report_csv,
          assembly_software_file = assembly_software_file
    }

    output {
        File sequencing_results_csv = results_table.sequencing_results_csv
        File wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
        File assembly_software_tsv = assembly_software_file
    }
}

task results_table {

    input {
      Array[String] sample_name
      File ev-d68_seq_results_py
      Array[File] cov_out
      Array[File] percent_cvg_csv
      String project_name
      File assembly_software_file
      File terra_data_table_path

    }

    command <<<
    python ~{ev-d68_seq_results_py} \
        --sample_name_array "~{write_lines(sample_name)}" \
        --cov_out_files "~{write_lines(cov_out)}" \
        --percent_cvg_files "~{write_lines(percent_cvg_csv)}" \
        --project_name "~{project_name}" \
        --assembly_software_file "~{assembly_software_file}" \
        --terra_data_table_path "~{terra_data_table_path}" \

    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"
        File wgs_horizon_report_csv = "~{project_name}_wgs_horizon_report.csv"
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
        File wgs_horizon_report_csv
        File assembly_software_file
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{wgs_horizon_report_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{assembly_software_file} ~{outdirpath}/summary_results/
    >>>

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
