# -----------------------------------------------------
# Virus abundance
# -----------------------------------------------------
import pandas as pd
import os


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
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus abundance rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 Align reads to virus catalog
# -----------------------------------------------------
# Align reads to virus catalog using bowtie2
rule build_viruses_bowtie2db:
    input:
        viruses=results
        + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
    output:
        results + "07_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog.1.bt2",
    params:
        db=results + "07_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input.viruses} {params.db} --threads {threads}
        """


# Align reads to virus catalog using bowtie2
rule align_reads_to_viruses:
    input:
        R1=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_2.fastq",
        R1S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_1.fastq",
        R2S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_2.fastq",
        db=results + "07_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog.1.bt2",
    output:
        results
        + "07_VIRUS_ABUNDANCE/01_align_viruses/bam_files/{group_assembly_sample}.bam",
    params:
        db=results + "07_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog",
        sam=results + "07_VIRUS_ABUNDANCE/01_align_viruses/{group_assembly_sample}.sam",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -U {input.R1S},{input.R2S} \
        -S {params.sam}

        # convert sam to bam
        samtools view -S -b {params.sam} > {output}
        rm {params.sam}
        """


rule prepare_read_counts_file:
    input:
        results + "01_READ_PREPROCESSING/read_preprocessing_report.csv",
    output:
        results + "07_VIRUS_ABUNDANCE/02_virus_abundance/read_counts.tsv",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/07_read_counts.py.ipynb"


# determine which viruses are present in the sample
rule metapop:
    input:
        bam=expand(
            results
            + "07_VIRUS_ABUNDANCE/01_align_viruses/bam_files/{group_assembly_sample}.bam",
            group_assembly_sample=group_assembly_sample,
        ),
        read_counts=results + "07_VIRUS_ABUNDANCE/02_virus_abundance/read_counts.tsv",
        viruses=results
        + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
    output:
        results + "07_VIRUS_ABUNDANCE/test",
    params:
        bam_dir=results + "07_VIRUS_ABUNDANCE/01_align_viruses/bam_files/",
        viruses_dir=results + "06_VIRUS_QUALITY/02_quality_filter/",
        out_dir=results + "07_VIRUS_ABUNDANCE/02_virus_abundance/",
        min_breadth=config["virus_abundance"]["min_breadth"],
        min_length=config["virus_abundance"]["min_length"],
        min_depth=config["virus_abundance"]["min_depth"],
    conda:
        "../envs/metapop.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        pip install metapop

        # run metapop to identify viruses present in samples
        metapop --input_samples {params.bam_dir} \
        --norm {input.read_counts} \
        --reference {params.viruses_dir} \
        --output {params.out_dir} \
        --min_cov {params.min_breadth} \
        --minimum_bases_for_detection {params.min_length} \
        --min_dep {params.min_depth} \
        --no_micro
        """
