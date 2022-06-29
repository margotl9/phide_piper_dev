# -----------------------------------------------------
# Virome enrichment (Only runs if data_type: "reads")
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")
assembly_sample = (
    samples_df["group"] + "_" + samples_df["assembly"] + "_" + samples_df["sample"]
)


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Enrichment rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 ViromeQC
# -----------------------------------------------------
# set up viromeqc
rule build_viromeqc:
    output:
        amphora=resources + "viromeqc/index/amphora_bacteria.dmnd",
        lsu=resources + "viromeqc/index/SILVA_132_LSURef_tax_silva.clean.1.bt2",
        ssu=resources + "viromeqc/index/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2",
    params:
        viromeqc_dir=resources + "viromeqc",
        index_dir=resources + "viromeqc/index",
        viromeqc_tmp=resources + "viromeqc/tmp/",
    shell:
        """
        # clone viromeqc
        rm -rf {params.viromeqc_dir}
        git clone --recurse-submodules https://github.com/SegataLab/viromeqc.git {params.viromeqc_dir}

        # install the viromeqc databases
        mkdir {params.index_dir}
        cd {params.index_dir}
        wget -O amphora_markers.zip https://zenodo.org/record/4020594/files/amphora_markers.zip?download=1
        wget -O SILVA_132_LSURef_tax_silva_clean.zip https://zenodo.org/record/4020594/files/SILVA_132_LSURef_tax_silva_clean.zip?download=1
        wget -O SILVA_132_SSURef_Nr99_tax_silva.clean.zip https://zenodo.org/record/4020594/files/SILVA_132_SSURef_Nr99_tax_silva.clean.zip?download=1

        # unzip databases
        unzip amphora_markers.zip
        unzip SILVA_132_LSURef_tax_silva_clean.zip
        unzip SILVA_132_SSURef_Nr99_tax_silva.clean.zip

        # make tmp dir
        mkdir {params.viromeqc_tmp}
        """


# determine virus enrichment with viromeqc
rule viromeqc:
    input:
        amphora=resources + "viromeqc/index/amphora_bacteria.dmnd",
        lsu=resources + "viromeqc/index/SILVA_132_LSURef_tax_silva.clean.1.bt2",
        ssu=resources + "viromeqc/index/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2",
        R1=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_2.fastq",
        R1S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_1.fastq",
        R2S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_2.fastq",
    output:
        results + "02_VIRUS_ENRICHMENT/01_viromeqc/{group_assembly_sample}/vqc.tsv",
    params:
        viromeqc_script=resources + "viromeqc/viromeQC.py",
        temp=resources + "viromeqc/tmp/{group_assembly_sample}/",
        extra_arguments=config["virus_enrichment"]["viromeqc_arguments"],
    log:
        results + "00_LOGS/02_virus_enrichment_{group_assembly_sample}.viromeqc.log",
    conda:
        "../envs/viromeqc.yml"
    threads: config["virus_enrichment"]["viromeqc_threads"]
    shell:
        """
        # make dir to act as tmp
        mkdir {params.temp}

        # determine virome enrichment using viromeqc
        {params.viromeqc_script} \
        --input {input.R1} {input.R2} {input.R1S} {input.R2S} \
        --output {output} \
        --bowtie2_threads {threads} \
        --diamond_threads {threads} \
        --tempdir {params.temp} \
        {params.extra_arguments} > {log} 2>&1

        # add sample column to each vqc output
        s={wildcards.group_assembly_sample}
        sed -i "s/$/\t$s/" {output}
        sample="sample"
        sed -i "1s/$s/$sample/" {output}

        # remove tmp dir
        rm -rf {params.temp}
        """


# combine viromeqc results for all samples
checkpoint combine_viromeqc_results_across_samples:
    input:
        expand(
            results
            + "02_VIRUS_ENRICHMENT/01_viromeqc/{group_assembly_sample}/vqc.tsv",
            group_assembly_sample=group_assembly_sample,
        ),
    output:
        results + "02_VIRUS_ENRICHMENT/virus_enrichment_report.tsv",
    shell:
        """
        # combine all viromeqc outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# Analyze viral enrichment
# -----------------------------------------------------
# analyze viromeqc results to visualize virus enrichment
rule virus_enrichment_analysis:
    input:
        results + "02_VIRUS_ENRICHMENT/virus_enrichment_report.tsv",
    output:
        report(
            results + "02_VIRUS_ENRICHMENT/virus_enrichment_figure.png",
            caption="../report/virus_enrichment_analysis.rst",
            category="Step 02: Virus enrichment",
        ),
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/02_virus_enrichment_analysis.py.ipynb"
