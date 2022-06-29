# -------------------------------------
# Read Preprocessing Module
# -------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")
group_assembly_sample = (
    samples_df["group"] + "_" + samples_df["assembly"] + "_" + samples_df["sample"]
)

# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -------------------------------------
# Preprocessing Rules
# -------------------------------------
# -----------------------------------------------------
# 00 Input
# -----------------------------------------------------
# symlink input paths to new paths
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[
            (
                samples_df["group"]
                + "_"
                + samples_df["assembly"]
                + "_"
                + samples_df["sample"]
                + "_"
                + samples_df["replicate"]
            )
            == wildcards.group_assembly_sample_replicate
        ]["R1"].iloc[0],
        R2=lambda wildcards: samples_df[
            (
                samples_df["group"]
                + "_"
                + samples_df["assembly"]
                + "_"
                + samples_df["sample"]
                + "_"
                + samples_df["replicate"]
            )
            == wildcards.group_assembly_sample_replicate
        ]["R2"].iloc[0],
    output:
        R1=results + "00_INPUT/{group_assembly_sample_replicate}_R1.fastq.gz",
        R2=results + "00_INPUT/{group_assembly_sample_replicate}_R2.fastq.gz",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# identify replicates
samples_df["group_assembly_sample"] = (
    samples_df["group"] + "_" + samples_df["assembly"] + "_" + samples_df["sample"]
)
sam_rep = samples_df[["group_assembly_sample", "replicate"]]
sam_rep_dict = sam_rep.set_index("group_assembly_sample").to_dict()["replicate"]


# -----------------------------------------------------
# 01 Merge replicates
# -----------------------------------------------------
# merge replicate files into single file
rule merge_replicates:
    input:
        R1=lambda wildcards: expand(
            results + "00_INPUT/{{group_assembly_sample}}_{replicate}_R1.fastq.gz",
            replicate=sam_rep_dict[wildcards.group_assembly_sample],
        ),
        R2=lambda wildcards: expand(
            results + "00_INPUT/{{group_assembly_sample}}_{replicate}_R2.fastq.gz",
            replicate=sam_rep_dict[wildcards.group_assembly_sample],
        ),
    output:
        R1=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_assembly_sample}_R1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_assembly_sample}_R2.fastq.gz",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 02 Gunzip files
# -----------------------------------------------------
# gunzip merged files
rule gunzip_merged_reads:
    input:
        R1=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_assembly_sample}_R1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_assembly_sample}_R2.fastq.gz",
    output:
        R1=results
        + "01_READ_PREPROCESSING/02_gunzip_reads/{group_assembly_sample}_R1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/02_gunzip_reads/{group_assembly_sample}_R2.fastq",
    shell:
        """
        gunzip -c {input.R1} > {output.R1}
        gunzip -c {input.R2} > {output.R2}
        """


# -----------------------------------------------------
# 03 Clumpify
# -----------------------------------------------------
# run clumpify to deduplicate reads
rule clumpify:
    input:
        R1=results
        + "01_READ_PREPROCESSING/02_gunzip_reads/{group_assembly_sample}_R1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/02_gunzip_reads/{group_assembly_sample}_R2.fastq",
    output:
        R1=results
        + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}_R1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}_R2.fastq",
        log=results + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}.log",
    params:
        extra_args=config["read_preprocessing"]["clumpify_args"],
    log:
        results + "00_LOGS/01_read_preprocessing_{group_assembly_sample}.clumpify.log",
    conda:
        "../envs/clumpify.yml"
    shell:
        """
        # run clumpify
        clumpify.sh \
        in={input.R1} \
        in2={input.R2} \
        out={output.R1} \
        out2={output.R2} \
        {params.extra_args} > {output.log} 2>&1

        cp {output.log} {log}
        """


# -----------------------------------------------------
# 04 KneadData
# -----------------------------------------------------
# build kneaddata bowtie2 database
rule download_kneaddata_database:
    output:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
    params:
        kneaddata_db=resources + "kneaddata/",
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # download human genome reference to desired directory
        kneaddata_database --download human_genome bowtie2 {params.kneaddata_db}
        """


# Quality filter and remove human reads with kneaddata
rule kneaddata:
    input:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=results
        + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}_R1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}_R2.fastq",
    output:
        log=results + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}.log",
        R1=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_2.fastq",
        R1S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_1.fastq",
        R2S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_2.fastq",
    params:
        output_dir=results + "01_READ_PREPROCESSING/04_kneaddata/",
        human_db=resources + "kneaddata/",
        extra_args=config["read_preprocessing"]["kneaddata_args"],
        prefix="{group_assembly_sample}",
    log:
        results + "00_LOGS/01_read_preprocessing_{group_assembly_sample}.kneaddata.log",
    conda:
        "../envs/kneaddata.yml"
    threads: config["read_preprocessing"]["kneaddata_threads"]
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.output_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        {params.extra_args}

        cp {output.log} {log}
        """


# -----------------------------------------------------
# Read preprocessing analysis
# -----------------------------------------------------
# Determine clumpify read counts
rule clumpify_read_counts:
    input:
        expand(
            results + "01_READ_PREPROCESSING/03_clumpify/{group_assembly_sample}.log",
            group_assembly_sample=group_assembly_sample,
        ),
    output:
        results + "01_READ_PREPROCESSING/03_clumpify/combined_read_counts.tsv",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/01_read_preprocessing_clumpify_read_counts.py.ipynb"


# determine read counts using kneaddata utils
rule kneaddata_read_counts:
    input:
        expand(
            results + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}.log",
            group_assembly_sample=group_assembly_sample,
        ),
    output:
        results + "01_READ_PREPROCESSING/04_kneaddata/combined_read_counts.tsv",
    params:
        log_dir=results + "01_READ_PREPROCESSING/04_kneaddata",
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        # generate read counts from kneaddata log files
        kneaddata_read_count_table \
        --input {params.log_dir} \
        --output {output}
        """


# Visualize read counts
rule read_count_analysis:
    input:
        clumpify=results + "01_READ_PREPROCESSING/03_clumpify/combined_read_counts.tsv",
        kneaddata=results
        + "01_READ_PREPROCESSING/04_kneaddata/combined_read_counts.tsv",
    output:
        figure=report(
            results + "01_READ_PREPROCESSING/read_preprocessing_figure.png",
            caption="../report/01_read_preprocessing_analysis.rst",
            category="Step 01: Read preprocessing",
        ),
        report=results + "01_READ_PREPROCESSING/read_preprocessing_report.csv",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/01_read_preprocessing_analysis.py.ipynb"
