version 1.0

workflow EVD68_illumina_pe_assembly {

    input {
        String    sample_name
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    evd68_genome
        File    evd68_gff
        String  project_name
        String out_dir

        # python scripts
        File    evd68_calc_percent_coverage_py
        File    evd68_concat_assembly_software_illumina_py  
    }

    # secret variables
    String outdirpath = sub(out_dir, "/$", "")

    call seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
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

    call create_software_assembly_file {
        input:
            evd68_concat_assembly_software_illumina_py = evd68_concat_assembly_software_illumina_py,
            project_name = project_name,
            bwa_version = align_reads.assembler_version,
            ivar_version = ivar_consensus.ivar_version
            
    }

    call transfer_outputs {
        input:
            fastqc_raw1_html = fastqc_raw1_html,
            fastqc_raw2_html = fastqc_raw2_html,
            fastqc_clean1_html = fastqc_clean1_html,
            fastqc_clean2_html = fastqc_clean2_html,
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
            out_dir = outdirpath
    }

    output {
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
        File assembly_software_file = create_software_assembly_file.assembly_software_file
        String bwa_version = align_reads.assembler_version
        String ivar_version = ivar_consensus.ivar_version
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

    command {

        seqyclean -minlen 70 -qual 20 20 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_name}_clean

    }

    output {

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"

    }

    runtime {
        cpu:    2
        memory:    "6 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/seqyclean:1.10.09"
    }
}

task fastqc {
    input {

        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command {

        fastqc --outdir $PWD ${fastq_1} ${fastq_2}

    }

    output {

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc2_html = "${fastq2_name}_fastqc.html"

    }

    runtime {
        cpu:    1
        memory:    "2 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/fastqc:0.11.9"
    }
}

task align_reads {

    input {

        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    command {

        echo bwa 0.7.17-r1188 > VERSION
        
        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_name}_aln.sorted.bam
        samtools index ./${sample_name}_aln.sorted.bam

    }

    output {

        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION")

    }

    runtime {
        cpu:    2
        memory:    "12 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "broadinstitute/viral-core:latest"
    }
}

task ivar_trim {

    input {

        File primers
        File bam
        String sample_name
    }

    command {

        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam

    }

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
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_var {

    input {

        String sample_name
        File ref
        File gff
        File bam
    }

    command {

        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}

    }

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
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {

    input {

        String sample_name
        File ref
        File bam
    }

    command <<<


        ivar version | awk '/version/ {print $3}' | tee VERSION

        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{ref} ~{bam} | \
        ivar consensus -p ~{sample_name}_consensus -q 20 -t 0.6 -m 10

    >>>

    output {

        File consensus_out = "${sample_name}_consensus.fa"
        String ivar_version = read_string("VERSION")

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "andersenlabapps/ivar:1.3.1"
    }
}

task bam_stats {

    input {

        String sample_name
        File bam
        File bai
    }

    command <<<

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

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.16"
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

task create_software_assembly_file {
    meta {
        description: "pull assembly software into a sinlge tsv file"
    }

    input {
        File evd68_concat_assembly_software_illumina_py
        String bwa_version
        String ivar_version
        String project_name
    }

    command <<<

        python ~{evd68_concat_assembly_software_illumina_py} \
        --project_name "~{project_name}" \
        --bwa_version "~{bwa_version}" \
        --ivar_version "~{ivar_version}"

    >>>

    output {
        File assembly_software_file = '~{project_name}_assembly_software.tsv'
    }

    runtime {

      docker: "mchether/py3-bio:v4"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"

    }
}

task transfer_outputs {
    input {
        String out_dir
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
    }

    String out_dir_path = sub('${out_dir}', "/$", "")

    command <<<
        
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