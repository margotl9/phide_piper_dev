# -----------------------------------------------------
# Read abundance
# -----------------------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv("config/samples.tsv", sep="\t")
assembly = samples_df["assembly"]
assembly_sample = (samples_df["assembly"] + "_" + samples_df["sample"])


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Read abundance rules
# -----------------------------------------------------
# symlink input paths to new paths (abundance and sample)
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[
            (samples_df["assembly"] + "_" + samples_df["sample"])
            == wildcards.assembly_sample
            ].iloc[0]["R1"],
        R2=lambda wildcards: samples_df[
            (samples_df["assembly"] + "_" + samples_df["sample"])
            == wildcards.assembly_sample
            ].iloc[0]["R2"],
    output:
        R1=results+"00_INPUT/{assembly_sample}_paired_1.fastq.gz",
        R2=results+"00_INPUT/{assembly_sample}_paired_2.fastq.gz",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """



# -----------------------------------------------------
# 01 Align reads to virus catalog using Bowtie2
# -----------------------------------------------------
# Align reads to virus catalog using bowtie2
rule build_viruses_bowtie2db:
    input:
        config["virus_db"]
    output:
        results+"READ_ABUNDANCE/01_bowtie2/virus_catalogs/virus_catalog.1.bt2",
    params:
        db=results+"READ_ABUNDANCE/01_bowtie2/virus_catalogs/virus_catalog", 
		extra_args=config["read_abundances"]["extra_args"]
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads {threads}
        {params.extra_args}
        """



# Align reads to virus catalog using bowtie2
rule bowtie2_align_reads_to_viruses:
    input:
        R1=results+"00_INPUT/{assembly_sample}_paired_1.fastq.gz",
        R2=results+"00_INPUT/{assembly_sample}_paired_2.fastq.gz",
        # R1S=results
        # + "01_READ_PREPROCESSING/03_kneaddata/{assembly_sample}_unmatched_1.fastq",
        # R2S=results
        # + "01_READ_PREPROCESSING/03_kneaddata/{assembly_sample}_unmatched_2.fastq",
        db=results+"READ_ABUNDANCE/01_bowtie2/virus_catalogs/virus_catalog.1.bt2",
    output:
        results
        + "READ_ABUNDANCE/01_bowtie2/bam_files/{assembly_sample}.bam",
    params:
        db=results+"READ_ABUNDANCE/01_bowtie2/virus_catalogs/virus_catalog",
        sam=results+"READ_ABUNDANCE/01_bowtie2/{assembly_sample}.sam",
        extra_args=config["read_abundances"]["extra_args"]

    log:
        results+"READ_ABUNDANCE/01_bowtie2/logs/{assembly_sample}_log",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    ##uncomment below and add to shell when using R1S and R2S 
    # -U {input.R1S},{input.R2S} \
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {params.sam} \
        > {log} 2>&1 \
        {params.extra_args}
 
        # convert sam to bam
        samtools view -S -b {params.sam} > {output} 
        rm {params.sam}
        """


# ----------------------------------------------------------
# 02 Align reads to virus database using Kraken2 and Bracken
# ----------------------------------------------------------
rule customize_virus_headers:
    input:
        genomes=config["virus_db"],
        metadata=config["virus_db_meta"],
    output:
        results
        +"READ_ABUNDANCE/02_kraken2_bracken/kraken_formatted_mgv.fasta",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/customize_virus_headers.py.ipynb"



# Align reads to virus catalog using kraken2
rule build_viruses_kraken2db:
    input:
        results
        +"READ_ABUNDANCE/02_kraken2_bracken/kraken_formatted_mgv.fasta",
    output:
        results
        +"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/hash.k2d",
        #config["read_abundances"]["kraken_db"] + "hash.k2d",
    params:
        db=results+"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/"
        #db=config["read_abundances"]["kraken_db"]
        extra_args=config["read_abundances"]["extra_args"]

    conda:
        "../envs/kraken2.yml"
    shell:
        """
        kraken2-build --download-taxonomy --db {params.db}
        kraken2-build --add-to-library {input} --db {params.db}
        kraken2-build --build --db {params.db}
        {params.extra_args}
        """



# Align reads to virus catalog using kraken2
rule kraken2_align_reads_to_viruses:
    input:
        db=results+"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/hash.k2d",
        #db=config["read_abundances"]["kraken_db"] + "hash.k2d",
        R1=results+"00_INPUT/{assembly_sample}_paired_1.fastq.gz",
        R2=results+"00_INPUT/{assembly_sample}_paired_2.fastq.gz",
    output:
        classification=results
        +"READ_ABUNDANCE/02_kraken2_bracken/kraken/{assembly_sample}_kraken2.kraken",
        report=results
        +"READ_ABUNDANCE/02_kraken2_bracken/reports/{assembly_sample}_kraken2.kreport",
    params:
        db=results+"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/",
        #db=config["read_abundances"]["kraken_db"],
        extra_args=config["read_abundances"]["extra_args"]

    conda:
        "../envs/kraken2.yml"
    shell:
        """
        kraken2 --paired {input.R1} {input.R2} \
        --db {params.db} \
        --report {output.report} > {output.classification} \
        {params.extra_args}
        """



# Build custom bracken database and align reads
rule bracken_build:
    input:
        results
        +"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/hash.k2d",
        #config["read_abundances"]["kraken_db"] + "hash.k2d",
    output:
        results
        + "READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/database150mers.kmer_distrib",
        #config["read_abundances"]["kraken_db"] + "database150mers.kmer_distrib",
    params:
        db=results
        + "READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/",
        #db=config["read_abundances"]["kraken_db"],
        extra_args=config["read_abundances"]["extra_args"]

    threads: config["virus_abundance"]["metapop_threads"]
    conda:
        "../envs/bracken.yml"
    shell:
        """
        bracken-build -d {params.db} -t {threads} -k 35 -l 150 \
        {params.extra_args}
        """



# Determine which viruses are present and their abundances
rule bracken:
    input:
        db=results
        + "READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/database150mers.kmer_distrib",
        #db=config["read_abundances"]["kraken_db"] + "database150mers.kmer_distrib",
        report=results+"READ_ABUNDANCE/02_kraken2_bracken/reports/{assembly_sample}_kraken2.kreport",
    output:
        # report=results+"READ_ABUNDANCE/02_kraken2_bracken/reports/{assembly_sample}_kraken2_bracken_families.kreport",
        abundance=results + "READ_ABUNDANCE/02_kraken2_bracken/bracken/{assembly_sample}_bracken_abundances.bracken",
    params:
        db=results+"READ_ABUNDANCE/02_kraken2_bracken/mgv_kraken2db/",
        #db=config["read_abundances"]["kraken_db"],
        level=config["virus_abundance"]["taxon_level"]
        extra_args=config["read_abundances"]["extra_args"]

    threads: config["virus_abundance"]["metapop_threads"]
    conda:
        "../envs/bracken.yml"
    shell:
        """
        bracken -d {params.db} \
        -i {input.report} \
        -o {output.abundance} \
        -r 150 \
        -l {params.level} \
        -t {threads} \
        {params.extra_args}
        """




# -----------------------------------------------------
# 03 Metapop
# -----------------------------------------------------
# rule prepare_read_counts_file:
#     input:
#         results + "01_READ_PREPROCESSING/read_preprocessing_report.csv",
#     output:
#         results + "READ_ABUNDANCE/03_metapop/read_counts.tsv",
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/07_read_counts.py.ipynb"

rule read_counts_from_bam:
    input:
        results
        + "READ_ABUNDANCE/01_bowtie2/bam_files/{assembly_sample}.bam",
    output:
        results+"READ_ABUNDANCE/03_metapop/{assembly_sample}_read_counts.tsv",
    conda:
        "../envs/kneaddata.yml"
    shell: 
        """
        echo -e {wildcards.assembly_sample}”\t”$(samtools view -c {input} > {output})
        """

rule combine_read_counts_across_samples:
    input:
        expand(results + "READ_ABUNDANCE/03_metapop/{assembly_sample}_read_counts.tsv", assembly_sample=assembly_sample),
    output:
        results + "READ_ABUNDANCE/03_metapop/combined_read_counts.tsv",
    conda:
        "../envs/kneaddata.yml"
    shell:
        """
        cat {input} > {output}
        """

# determine which viruses are present in the sample
rule metapop:
    input:
        bam=expand(
            results
            +"READ_ABUNDANCE/01_bowtie2/bam_files/{assembly_sample}.bam",
            assembly_sample=assembly_sample),
        read_counts=results + "READ_ABUNDANCE/03_metapop/combined_read_counts.tsv",
        viruses=config["virus_db"]
    output:
        results+"READ_ABUNDANCE/test",
    params:
        bam_dir=results + "READ_ABUNDANCE/01_bowtie2/bam_files/",
        viruses_dir=resources + "mgv_db/", #config["read_abundances"]["virus_dir"],
        out_dir=results + "READ_ABUNDANCE/03_metapop/",
        min_breadth=config["virus_abundance"]["min_breadth"],
        min_length=config["virus_abundance"]["min_length"],
        min_depth=config["virus_abundance"]["min_depth"],
        extra_args=config["read_abundances"]["extra_args"],

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
        --no_micro \
        --threads {threads} \
        {params.extra_args}
        """



# -----------------------------------------------------
# 04 InStrain profile
# -----------------------------------------------------
rule contig_from_virus_db:
    input:
        config["virus_db"]
    output:
        resources+"mgv_db/contig_mgv.fna"
        #config["read_abundances"]["virus_dir"] + "contig_mgv.fna"
    conda:
        "../envs/instrain.yml"
    shell:
        """
        prodigal -i {input} \
        -o {output}
        """

rule sort_bam:
    input: 
        bam=results+"READ_ABUNDANCE/01_bowtie2/bam_files/{assembly_sample}.bam",
    output:
        sort = results+"READ_ABUNDANCE/01_bowtie2/sorted_bam_files/{assembly_sample}_sorted.bam",
    conda:
        "../envs/instrain.yml"
    shell:
        """
        samtools sort {input.bam} -o {output.sort}
        """

rule inStrain_profile:
    input:
        # fna=resources+"mgv_db/contig_mgv.fna",
        #config["read_abundances"]["virus_dir"] + "contig_mgv.fna"
        fasta=config["virus_db"],
        bam=results+"READ_ABUNDANCE/01_bowtie2/sorted_bam_files/{assembly_sample}_sorted.bam",
    output:
        results+"READ_ABUNDANCE/04_instrain_profile/{assembly_sample}_profile.IS/output/{assembly_sample}_profile.IS_linkage.tsv",
    params:
        out_dir=results+"READ_ABUNDANCE/04_instrain_profile/{assembly_sample}_profile.IS",
        min_breadth=config["virus_abundance"]["min_breadth"],
        extra_args=config["read_abundances"]["extra_args"]

    conda:
        "../envs/instrain.yml"
    # -g {input.fna} \
    shell:
        """
        inStrain profile \
        {input.bam} \
        {input.fasta} \
        -o {params.out_dir} \
        --skip_mm_profiling \
        --skip_genome_wide \
        --min_cov {params.min_breadth} \
        {params.extra_args}
        """



# -----------------------------------------------------
# 04 InStrain compare
# -----------------------------------------------------
rule inStrain_compare:
    input:
        expand(results+"READ_ABUNDANCE/04_instrain_profile/{assembly_sample}_profile.IS/output/{assembly_sample}_profile.IS_linkage.tsv", assembly_sample=assembly_sample)
    output:
        results+"READ_ABUNDANCE/05_instrain_compare/compare.IS/output/compare.IS_comparisonsTable.tsv",
    params:
        bam=expand(results+"READ_ABUNDANCE/01_bowtie2/sorted_bam_files/{assembly_sample}_sorted.bam", assembly_sample=assembly_sample),
        profile=expand(results+"READ_ABUNDANCE/04_instrain_profile/{assembly_sample}_profile.IS", assembly_sample=assembly_sample),
        out_dir=results+"READ_ABUNDANCE/05_instrain_compare/compare.IS",
        extra_args=config["read_abundances"]["extra_args"]

    conda:
        "../envs/instrain.yml"
    shell:
        """
        inStrain compare \
        -i {params.profile} \
        -o {params.out_dir} \
        -bams {params.bam} \
        {params.extra_args}
        """


