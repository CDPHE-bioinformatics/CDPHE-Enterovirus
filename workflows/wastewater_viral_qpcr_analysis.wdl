version 1.0

# import workflow version capture task
import "../tasks/hostile_task.wdl" as hostile_task

workflow wastewater_viral_qpcr_analysis {

    input {
        String    sample_name
        File    fastq_1
        File    fastq_2
        File    primer_bed
        File    adapters_and_contaminants
        File    ref_genome
        Boolean    scrub_reads
        Array[File]?    scrub_genome_index
        String    project_name
        String  qpcr_region
        String out_dir

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
            ref = ref_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_trim {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call bam_stats {
        input:
            sample_name = sample_name,
            bam = ivar_trim.trimsort_bam,
            bai = ivar_trim.trimsort_bamindex,
            qpcr_region = qpcr_region
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
            sorted_bam = sorted_bam,
            sorted_bamindex = sorted_bamindex,
            trimsort_bam = trimsort_bam,
            trimsort_bamindex = trimsort_bamindex,
            cov_out = cov_out,
            cov_out_10x = cov_out_10x,
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
        File sorted_bam = align_reads.sorted_bam
        File sorted_bamindex = align_reads.sorted_bamindex
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex
        File cov_out = bam_stats.cov_out
        File cov_out_10x = bam_stats.cov_out_10x
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

        seqyclean -minlen 25 -qual 30 30 -gz -1 ~{fastq_1} -2 ~{fastq_2} -c ~{contam} -o ~{sample_name}_clean

    >>>

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

    String docker = "staphb/fastqc:0.12.1"

    command <<<

        fastqc --outdir $PWD ~{fastq_1} ~{fastq_2}

    >>>

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

        bwa index -p reference.fasta -a is ~{ref}
        bwa mem -t 2 reference.fasta ~{fastq_1} ~{fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./~{sample_name}_aln.sorted.bam
        samtools index ./~{sample_name}_aln.sorted.bam

    >>>

    output {

        File sorted_bam = "${sample_name}_aln.sorted.bam"
        File sorted_bamindex = "${sample_name}_aln.sorted.bam.bai"

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

        ivar trim -e -i ~{bam} -b ~{primers} -p ~{sample_name}_trim.bam
        samtools sort ~{sample_name}_trim.bam -o ~{sample_name}_trim.sort.bam
        samtools index ~{sample_name}_trim.sort.bam

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

task bam_stats {

    input {

        String sample_name
        String qpcr_region
        File bam
        File bai
    }

    String docker = "staphb/samtools:1.23"

    command <<<
    
        samtools coverage -r ~{qpcr_region} -o ~{sample_name}_coverage.txt ~{bam}
        samtools coverage --min-depth 10-r ~{qpcr_region} -o ~{sample_name}_coverage_10x.txt ~{bam}

    >>>

    output {

        File cov_out  = "${sample_name}_coverage.txt"
        File cov_out_10x  = "${sample_name}_coverage_10x.txt"

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
        File sorted_bam
        File sorted_bamindex
        File trimsort_bam
        File trimsort_bamindex
        File cov_out
        File cov_out_10x
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
        gsutil -m cp ~{sorted_bam} ~{out_dir_path}/alignments/
        gsutil -m cp ~{sorted_bamindex} ~{out_dir_path}/alignments/
        gsutil -m cp ~{trimsort_bam} ~{out_dir_path}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{out_dir_path}/alignments/
        gsutil -m cp ~{cov_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{cov_out_10x} ~{out_dir_path}/bam_stats/

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