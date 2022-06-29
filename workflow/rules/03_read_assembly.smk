# -------------------------------------
# Read assembly (Only runs if data_type: "reads")
# -------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -------------------------------------
# Read assembly rules
# -------------------------------------
# -----------------------------------------------------
# 01 combine reads for coassembly
# -----------------------------------------------------
# identify replicates
samples_df["group_assembly"] = samples_df["group"] + "_" + samples_df["assembly"]
group_assemb = samples_df[["group_assembly", "sample"]]
group_assem_dict = group_assemb.set_index("group_assembly").to_dict()["sample"]


# combine enriched reads for assembly (or coassembly if specified)
rule merge_reads_for_assembly:
    input:
        R1=lambda wildcards: expand(
            results
            + "01_READ_PREPROCESSING/04_kneaddata/{{group_assembly}}_{sample}_paired_1.fastq",
            sample=group_assem_dict[wildcards.group_assembly],
        ),
        R2=lambda wildcards: expand(
            results
            + "01_READ_PREPROCESSING/04_kneaddata/{{group_assembly}}_{sample}_paired_2.fastq",
            sample=group_assem_dict[wildcards.group_assembly],
        ),
        R1S=lambda wildcards: expand(
            results
            + "01_READ_PREPROCESSING/04_kneaddata/{{group_assembly}}_{sample}_unmatched_1.fastq",
            sample=group_assem_dict[wildcards.group_assembly],
        ),
        R2S=lambda wildcards: expand(
            results
            + "01_READ_PREPROCESSING/04_kneaddata/{{group_assembly}}_{sample}_unmatched_2.fastq",
            sample=group_assem_dict[wildcards.group_assembly],
        ),
    output:
        R1=results + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_paired_1.fastq",
        R2=results + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_paired_2.fastq",
        R1S=results
        + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_unmatched_1.fastq",
        R2S=results
        + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_unmatched_2.fastq",
    shell:
        """
        # combine reads for coassembly
        cat {input.R1} > {output.R1}
        cat {input.R2} > {output.R2}
        cat {input.R1S} > {output.R1S}
        cat {input.R2S} > {output.R2S}
        """


# -----------------------------------------------------
# 02 metaSPAdes
# -----------------------------------------------------
# assemble reads using metaspades
rule metaspades:
    input:
        R1=results + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_paired_1.fastq",
        R2=results + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_paired_2.fastq",
        R1S=results
        + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_unmatched_1.fastq",
        R2S=results
        + "03_READ_ASSEMBLY/01_combine_reads/{group_assembly}_unmatched_2.fastq",
    output:
        results
        + "03_READ_ASSEMBLY/02_metaspades/{group_assembly}/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/02_metaspades/{group_assembly}",
        extra_args=config["read_assembly"]["metaspades_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{group_assembly}.metaspades.log",
    conda:
        "../envs/metaspades.yml"
    threads: config["read_assembly"]["metaspades_threads"]
    shell:
        """
        # assemble reads using metaspades
        spades.py \
        --meta \
        --pe1-1 {input.R1} \
        --pe1-2 {input.R2} \
        --pe1-s {input.R1S} \
        --pe1-s {input.R2S} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/spades.log {log}
        """


# -----------------------------------------------------
# 03 QUAST
# -----------------------------------------------------
# run quast to determine the quality of the assemblies
rule quast:
    input:
        results
        + "03_READ_ASSEMBLY/02_metaspades/{group_assembly}/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    output:
        results + "03_READ_ASSEMBLY/03_quast/{group_assembly}/transposed_report.tsv",
    params:
        output_dir=results + "03_READ_ASSEMBLY/03_quast/{group_assembly}",
        min_len=config["read_assembly"]["min_contig_length"],
        labels="{group_assembly}",
        extra_args=config["read_assembly"]["quast_arguments"],
    log:
        results + "00_LOGS/03_read_assembly_{group_assembly}.quast.log",
    conda:
        "../envs/quast.yml"
    shell:
        """
        # assembly analysis using quast
        metaquast.py \
        {input} \
        -o {params.output_dir} \
        --threads {threads} \
        --min-contig {params.min_len} \
        --contig-thresholds 0,1000,5000,10000,{params.min_len} \
        --labels {params.labels} \
        {params.extra_args}

        # copy spades.log to log file
        cp {params.output_dir}/quast.log {log}
        """


# -----------------------------------------------------
# 04 Contig length filter
# -----------------------------------------------------
# filter contigs based on contig length
rule contig_length_filter:
    input:
        results
        + "03_READ_ASSEMBLY/02_metaspades/{group_assembly}/"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    output:
        results
        + "03_READ_ASSEMBLY/04_contig_length_filter/{group_assembly}_"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/03_read_assembly_contig_length_filter.py.ipynb"


# Determine inputs for combine_quast_outputs
def get_enriched_assemblies(wildcards):
    vqc = pd.read_csv(
        checkpoints.combine_viromeqc_results_across_samples.get(**wildcards).output[0],
        sep="\t",
    )
    vqc_hq = vqc[
        vqc["total enrichmnet score"] > config["virus_enrichment"]["min_enrichment"]
    ]
    vqc_hq["group_assembly"] = vqc_hq["sample"].str.rpartition("_")[0]
    return expand(
        results + "03_READ_ASSEMBLY/03_quast/{group_assembly}/transposed_report.tsv",
        group_assembly=list(set(vqc_hq["group_assembly"])),
    )


# combine quast outputs
checkpoint combine_quast_outputs_across_samples:
    input:
        get_enriched_assemblies,
    output:
        results + "03_READ_ASSEMBLY/read_assembly_report.tsv",
    shell:
        """
        # combine quast reports for all assemblies, only keeping the header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# Analyze assemblies
# -----------------------------------------------------
# analyze quast results to visualize assembly quality
rule read_assembly_analysis:
    input:
        results + "03_READ_ASSEMBLY/read_assembly_report.tsv",
    output:
        report(
            results + "03_READ_ASSEMBLY/read_assembly_figure.png",
            caption="../report/read_assembly_analysis_contig_count.rst",
            category="Step 02: Read assembly",
        ),
    params:
        min_len=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/03_read_assembly_analysis.py.ipynb"
