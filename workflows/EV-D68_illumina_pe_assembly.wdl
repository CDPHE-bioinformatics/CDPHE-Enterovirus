version 1.0

# import workflow version capture task
import "../tasks/version_capture_task.wdl" as version_capture
import "../tasks/hostile_task.wdl" as hostile_task

workflow EVD68_illumina_pe_assembly {

    input {
        String    sample_name
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    evd68_genome
        File    evd68_gff
        Boolean scrub_reads
        Array[File]? scrub_genome_index
        String  project_name
        String out_dir

        # python scripts
        File    evd68_calc_percent_coverage_py
        File    version_capture_py
    }

    # secret variables
    String outdirpath = sub(out_dir, "/$", "")

    if (scrub_reads) {
      call hostile_task.hostile as hostile {
          input:
              fastq1 = fastq_1,
              fastq2 = fastq_2,
              genome_index = select_first([scrub_genome_index]),
              seq_method = "ILLUMINA"
      }
  }

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
            fastq_1 = select_first([hostile.fastq1_scrubbed, fastq_1]),
            fastq_2 = select_first([hostile.fastq2_scrubbed, fastq_2])
    }

    call fastqc as fastqc_raw {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call fastqc as fastqc_cleaned {
        input:
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call align_reads {
        input:
            sample_name = sample_name,
            ref = evd68_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_trim {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call ivar_var {
        input:
            sample_name = sample_name,
            ref = evd68_genome,
            gff = evd68_gff,
            bam = ivar_trim.trimsort_bam
    }

    call ivar_consensus {
        input:
            sample_name = sample_name,
            ref = evd68_genome,
            bam = ivar_trim.trimsort_bam
    }

    call bam_stats {
        input:
            sample_name = sample_name,
            bam = ivar_trim.trimsort_bam,
            bai = ivar_trim.trimsort_bamindex,
    }

    call rename_fasta {
        input:
            sample_name = sample_name,
            fasta = ivar_consensus.consensus_out
    }

    call calc_percent_cvg {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            evd68_calc_percent_coverage_py = evd68_calc_percent_coverage_py
    }

    call version_capture.workflow_version_capture  as workflow_version_capture{
        input:
    }

    Array[VersionInfo] version_array = [
        seqyclean.seqyclean_version_info,
        fastqc_cleaned.fastqc_version_info,
        align_reads.bwa_version_info,
        align_reads.samtools_version_info,
        ivar_consensus.ivar_version_info,
        ivar_consensus.samtools_version_info,
        bam_stats.samtools_version_info
    ]
    if (scrub_reads) {
        Array[VersionInfo] version_array_with_hostile = flatten([version_array, select_all([hostile.hostile_version_info])])
    }

    call version_capture.task_version_capture as task_version_capture {
        input:
            version_array = select_first([version_array_with_hostile, version_array]),
            workflow_name = "EV-D68_illumina_pe_assembly",
            workflow_version = workflow_version_capture.workflow_version,
            project_name = project_name,
            analysis_date = workflow_version_capture.analysis_date,
            version_capture_py = version_capture_py

    }

    call transfer_outputs {
        input:
            fastq1_scrubbed = hostile.fastq1_scrubbed,
            fastq2_scrubbed = hostile.fastq2_scrubbed,
            fastqc_raw1_html = fastqc_raw.fastqc1_html,
            fastqc_raw2_html = fastqc_raw.fastqc2_html,
            fastqc_clean1_html = fastqc_cleaned.fastqc1_html,
            fastqc_clean2_html = fastqc_cleaned.fastqc2_html,
            seqyclean_summary = seqyclean_summary,
            filtered_reads_1 = filtered_reads_1,
            filtered_reads_2 = filtered_reads_2,
            trimsort_bam = trimsort_bam,
            trimsort_bamindex = trimsort_bamindex,
            consensus = consensus,
            variants = variants,
            cov_out = cov_out,
            covhist_out = covhist_out,
            flagstat_out = flagstat_out,
            stats_out = stats_out,
            depth_out = depth_out,
            renamed_consensus = renamed_consensus,
            version_capture_illumina_pe_assembly = task_version_capture.version_capture_file,
            out_dir = outdirpath
    }

    output {
        Int? human_reads_removed = hostile.human_reads_removed
        Float? human_reads_removed_proportion = hostile.human_reads_removed_proportion
        File? fastq1_scrubbed = hostile.fastq1_scrubbed
        File? fastq2_scrubbed = hostile.fastq2_scrubbed
        File filtered_reads_1 = seqyclean.cleaned_1
        File filtered_reads_2 = seqyclean.cleaned_2
        File seqyclean_summary = seqyclean.seqyclean_summary
        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_clean1_html = fastqc_cleaned.fastqc1_html
        File fastqc_clean2_html = fastqc_cleaned.fastqc2_html
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex
        File variants = ivar_var.var_out
        File consensus = ivar_consensus.consensus_out
        File flagstat_out = bam_stats.flagstat_out
        File stats_out = bam_stats.stats_out
        File depth_out = bam_stats.depth_out
        File covhist_out = bam_stats.covhist_out
        File cov_out = bam_stats.cov_out
        File renamed_consensus = rename_fasta.renamed_consensus
        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv
        File version_capture_illumina_pe_assembly = task_version_capture.version_capture_file
        String transfer_date = transfer_outputs.transfer_date
    }
}

task seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/seqyclean:1.10.09"

    command <<<

        seqyclean -minlen 70 -qual 20 20 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_name}_clean

        # grab seqyclean version
        seqyclean -h | awk '/Version/ {print $2}' | tee VERSION

    >>>

    output {

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"

        VersionInfo seqyclean_version_info = object {
            software: "seqyclean",
            docker: docker,
            version: read_string("VERSION")
        }

    }

    runtime {
        cpu:    2
        memory:    "6 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task fastqc {
    input {

        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    String docker = "staphb/fastqc:0.11.9"

    command <<<

        fastqc --outdir $PWD ${fastq_1} ${fastq_2}

        # grab version
        fastqc --version | awk '/FastQC/ {print $2}' | tee VERSION

    >>>

    output {

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc2_html = "${fastq2_name}_fastqc.html"

        VersionInfo fastqc_version_info = object {
            software: "fastqc",
            docker: docker,
            version: read_string("VERSION")
        }

    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task align_reads {

    input {

        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    String docker = "quay.io/broadinstitute/viral-core:2.2.3"

    command <<<

        # echo bwa 0.7.17-r1188 > VERSION
        # grab version bwa and samtools versions
        bwa 2>&1 | awk '/Version/{print $2}' | tee VERSION_bwa
        samtools --version | awk '/samtools/ {print $2}' |tee VERSION_samtools

        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_name}_aln.sorted.bam
        samtools index ./${sample_name}_aln.sorted.bam

    >>>

    output {

        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"

        VersionInfo bwa_version_info = object {
            software: "bwa",
            docker: docker,
            version: read_string("VERSION_bwa")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_samtools")
        }

    }

    runtime {
        cpu:    2
        memory:    "12 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam

    >>>

    output {

        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_var {

    input {

        String sample_name
        File ref
        File gff
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}

    >>>

    output {

        File var_out = "${sample_name}_variants.tsv"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task ivar_consensus {

    input {

        String sample_name
        File ref
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<

        # grab ivar and samtools versions
        ivar version | awk '/version/ {print $3}' | tee VERSION_ivar
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_samtools

        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{ref} ~{bam} | \
        ivar consensus -p ~{sample_name}_consensus -q 20 -t 0.6 -m 10

    >>>

    output {

        File consensus_out = "${sample_name}_consensus.fa"

        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_ivar")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string ("VERSION_samtools")
        }

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task bam_stats {

    input {

        String sample_name
        File bam
        File bai
    }

    String docker = "staphb/samtools:1.16"

    command <<<

        # grab version
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION

        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}
        samtools depth -a -o ~{sample_name}_depth.txt ~{bam}

    >>>

    output {

        File flagstat_out  = "${sample_name}_flagstat.txt"
        File stats_out  = "${sample_name}_stats.txt"
        File covhist_out  = "${sample_name}_coverage_hist.txt"
        File cov_out  = "${sample_name}_coverage.txt"
        File depth_out  = "${sample_name}_depth.txt"

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    docker
    }
}

task rename_fasta {

    input {

        String sample_name
        File fasta
    }

    command <<<

        sed 's/>.*/>CO-CDPHE-~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa

    >>>

    output {

        File renamed_consensus  = "${sample_name}_consensus_renamed.fa"

    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task calc_percent_cvg {

    input {

        File fasta
        String sample_name
        File evd68_calc_percent_coverage_py

    }

    command {
        python ~{evd68_calc_percent_coverage_py} \
          --sample_name ~{sample_name} \
          --fasta_file ~{fasta}
      }
    output {

      File percent_cvg_csv  = "${sample_name}_consensus_cvg_stats.csv"

    }

    runtime {

      docker: "mchether/py3-bio:v1"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }

}

task transfer_outputs {
    input {
        String out_dir
        File? fastq1_scrubbed
        File? fastq2_scrubbed
        File fastqc_raw1_html
        File fastqc_raw2_html
        File fastqc_clean1_html
        File fastqc_clean2_html
        File seqyclean_summary
        File filtered_reads_1
        File filtered_reads_2
        File trimsort_bam
        File trimsort_bamindex
        File consensus
        File variants
        File cov_out
        File covhist_out
        File flagstat_out
        File stats_out
        File depth_out
        File renamed_consensus
        File version_capture_illumina_pe_assembly
    }

    String out_dir_path = sub('${out_dir}', "/$", "")

    command <<<

        gsutil -m cp ~{fastq1_scrubbed} ~{out_dir_path}/scrubbed_fastq/
        gsutil -m cp ~{fastq2_scrubbed} ~{out_dir_path}/scrubbed_fastq/
        gsutil -m cp ~{fastqc_raw1_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw2_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_clean1_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_clean2_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{seqyclean_summary} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{filtered_reads_1} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{filtered_reads_2} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{trimsort_bam} ~{out_dir_path}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{out_dir_path}/alignments/
        gsutil -m cp ~{consensus} ~{out_dir_path}/assemblies/
        gsutil -m cp ~{variants} ~{out_dir_path}/variants/
        gsutil -m cp ~{cov_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{flagstat_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{stats_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{depth_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{renamed_consensus} ~{out_dir_path}/assemblies/
        gsutil -m cp ~{version_capture_illumina_pe_assembly} ~{out_dir_path}/summary_results/

        transferdate=`date`
        echo $transferdate | tee TRANSFERDATE
    >>>

    output {
        String transfer_date = read_string("TRANSFERDATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
