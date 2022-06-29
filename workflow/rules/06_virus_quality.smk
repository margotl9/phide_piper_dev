# -----------------------------------------------------
# Virus quality (will always run)
# -----------------------------------------------------
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


# -----------------------------------------------------
# Virus quality rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 CheckV
# -----------------------------------------------------
# download checkv database
rule build_checkv:
    output:
        resources + "checkv/checkv-db-v1.2/README.txt",
    params:
        checkv_dir=resources + "checkv/",
    conda:
        "../envs/checkv.yml"
    shell:
        """
        # download checkv database
        checkv download_database {params.checkv_dir}
        """


# run checkv on viral contigs to determine genome quality
rule checkv:
    input:
        checkv_db=resources + "checkv/checkv-db-v1.2/README.txt",
        virus_contigs=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/cluster_representatives.fna",
    output:
        checkv_results=results + "06_VIRUS_QUALITY/01_checkv/quality_summary.tsv",
        checkv_proviruses=results + "06_VIRUS_QUALITY/01_checkv/proviruses.fna",
        checkv_viruses=results + "06_VIRUS_QUALITY/01_checkv/viruses.fna",
    params:
        checkv_dir=results + "06_VIRUS_QUALITY/01_checkv",
        checkv_db=resources + "checkv/checkv-db-v1.2",
    log:
        results + "00_LOGS/06_virus_quality.checkv.log",
    conda:
        "../envs/checkv.yml"
    threads: config["virus_quality"]["checkv_threads"]
    shell:
        """
        # run checkv to determine virus quality
        checkv end_to_end {input.virus_contigs} {params.checkv_dir} \
        -d {params.checkv_db} \
        -t {threads} > {log} 2>&1
        """


# -----------------------------------------------------
# 02 Quality filter viruses
# -----------------------------------------------------
rule quality_filter_viruses:
    input:
        checkv_proviruses=results + "06_VIRUS_QUALITY/01_checkv/proviruses.fna",
        checkv_viruses=results + "06_VIRUS_QUALITY/01_checkv/viruses.fna",
        checkv_results=results + "06_VIRUS_QUALITY/01_checkv/quality_summary.tsv",
    output:
        results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
    params:
        min_completeness=config["virus_quality"]["min_completeness"],
        min_viral_genes=config["virus_quality"]["min_viral_genes"],
        max_bacterial_genes=config["virus_quality"]["max_bacterial_genes"]
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/06_quality_filter_viruses.py.ipynb"


# -----------------------------------------------------
# Analyze virus quality
# -----------------------------------------------------
# analyze checkv results to visualize genome qualities
rule virus_quality_analysis:
    input:
        results + "06_VIRUS_QUALITY/01_checkv/quality_summary.tsv",
    output:
        report(
            results + "06_VIRUS_QUALITY/virus_quality_figure.png",
            caption="../report/06_virus_quality_analysis.rst",
            category="Step 06: Virus quality",
        ),
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/06_virus_quality_analysis.py.ipynb"
