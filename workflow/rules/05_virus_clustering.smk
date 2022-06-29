# -----------------------------------------------------
# Virus clustering
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
# Virus clustering rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 Combine viruses
# -----------------------------------------------------
# Determine inputs for combine_quast_outputs
def get_viral_contigs(wildcards):
    virus_identification_report = pd.read_csv(
        checkpoints.combine_reports_across_samples.get(**wildcards).output[0],
    )
    return expand(
        results
        + "04_VIRUS_IDENTIFICATION/07_combine_outputs/{group_assembly}/combined_viral_contigs.fasta",
        group_assembly=list(set(virus_identification_report["assembly"])),
    )


rule combine_viral_contigs:
    input:
        get_viral_contigs,
    output:
        results + "05_VIRUS_CLUSTERING/01_combine_viruses/combined_viruses.fna",
    shell:
        """
        cat {input} > {output}
        """


# -----------------------------------------------------
# 02 Dereplicate viruses
# -----------------------------------------------------
# dereplicate viral genomes
rule dereplicate_viruses:
    input:
        results + "05_VIRUS_CLUSTERING/01_combine_viruses/combined_viruses.fna",
    output:
        results + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/virus_replicates.tsv",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        blastdb=results + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/virus_blastdb",
        blast_tsv=results
        + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/viruses_blast.tsv",
        ani_tsv=results + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/viruses_ani.tsv",
        min_ani=config["virus_clustering"]["derep_min_ani"],
        min_tcov=config["virus_clustering"]["derep_min_tcov"],
        min_qcov=config["virus_clustering"]["derep_min_qcov"],
    conda:
        "../envs/mgv.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        # make a blast db from phage contigs
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {input} -db {params.blastdb} -out {params.blast_tsv} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90

        # calculate ani and af from blast results
        python {params.blastani_script} -i {params.blast_tsv} -o {params.ani_tsv}

        # cluster phage genomes based on 95% ani and 85% af
        python {params.cluster_script} --fna {input} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


# extract cluster representatives from clusters
rule extract_dereplicated_viruses:
    input:
        viruses=results + "05_VIRUS_CLUSTERING/01_combine_viruses/combined_viruses.fna",
        clusters=results
        + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/virus_replicates.tsv",
    output:
        results + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/dereplicated_viruses.fna",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/05_extract_cluster_representatives.py.ipynb"


# -----------------------------------------------------
# 03 Cluster viruses
# -----------------------------------------------------
# cluster viral genomes
rule cluster_viruses:
    input:
        results + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/dereplicated_viruses.fna",
    output:
        results + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_clusters.tsv",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        blastdb=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_blastdb",
        blast_tsv=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/viruses_blast.tsv",
        ani_tsv=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/viruses_ani.tsv",
        min_ani=config["virus_clustering"]["min_ani"],
        min_tcov=config["virus_clustering"]["min_tcov"],
        min_qcov=config["virus_clustering"]["min_qcov"],
    conda:
        "../envs/mgv.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        # make a blast db from phage contigs
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {input} -db {params.blastdb} -out {params.blast_tsv} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90

        # calculate ani and af from blast results
        python {params.blastani_script} -i {params.blast_tsv} -o {params.ani_tsv}

        # cluster phage genomes based on 95% ani and 85% af
        python {params.cluster_script} --fna {input} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


# align reads to virus db
rule align_reads_to_virusdb:
    input:
        R1=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_paired_2.fastq",
        R1S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_1.fastq",
        R2S=results
        + "01_READ_PREPROCESSING/04_kneaddata/{group_assembly_sample}_unmatched_2.fastq",
        virusdb=config["virus_db"],
    output:
        results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/bam_files/{group_assembly_sample}_virus_v_virusdb.bam",
    params:
        db=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/virusdb",
        sam=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/bam_files/{group_assembly_sample}_virus_v_virusdb.sam",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input.virusdb} {params.db} --threads {threads}

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
        rm -rf {params.sam}
        """


# build metapop
rule build_metapop:
    output:
        results + "05_VIRUS_CLUSTERING/03_cluster_viruses/metapop_build_complete",
    params:
        bam_dir=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/bam_files/",
    conda:
        "../envs/metapop.yml"
    shell:
        """
        # build metapop
        pip install metapop

        touch {output}
        """


# determine which viruses are present in the sample
rule identify_present_viruses:
    input:
        metapop=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/metapop_build_complete",
        bam=expand(
            results
            + "05_VIRUS_CLUSTERING/03_cluster_viruses/bam_files/{group_assembly_sample}_virus_v_virusdb.bam",
            group_assembly_sample=group_assembly_sample,
        ),
        virusdb=config["virus_db"],
    output:
        results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/MetaPop/03.Breadth_and_Depth/group1_assembly1_enriched1_virus_v_virusdb_breadth_and_depth.tsv",
    params:
        bam_dir=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/bam_files/",
        out_dir=results + "05_VIRUS_CLUSTERING/03_cluster_viruses/",
    conda:
        "../envs/metapop.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        # run metapop to identify viruses present in samples
        metapop --input_samples {params.bam_dir} \
        --reference {input.virusdb} \
        --genes /home/carsonjm/resources/virusdb/proteins/all_genomes_genes.fasta \
        --output {params.out_dir} \
        --no_micro --no_macro
        """


# filter to keep only blast hits
rule extract_virusdb_hits:
    input:
        virusdb=config["virus_db"],
        metapop=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/MetaPop/03.Breadth_and_Depth/group1_assembly1_enriched1_virus_v_virusdb_breadth_and_depth.tsv",
    output:
        results + "05_VIRUS_CLUSTERING/03_cluster_viruses/virusdb_blast_hits.fna",
    params:
        min_breadth=config["virus_clustering"]["min_virusdb_breadth"],
        min_length=config["virus_clustering"]["min_virusdb_length"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/05_extract_virusdb_hits.py.ipynb"


# cluster viral genomes with virusdb
rule cluster_viruses_with_virusdb:
    input:
        viruses=results
        + "05_VIRUS_CLUSTERING/02_dereplicate_viruses/dereplicated_viruses.fna",
        virusdb=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virusdb_blast_hits.fna",
    output:
        clusters=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_clusters.tsv",
        combined=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined.fna",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        blastdb=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_blastdb",
        blast_tsv=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_blast.tsv",
        ani_tsv=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_ani.tsv",
        min_ani=config["virus_clustering"]["min_ani"],
        min_tcov=config["virus_clustering"]["min_tcov"],
        min_qcov=config["virus_clustering"]["min_qcov"],
    conda:
        "../envs/mgv.yml"
    threads: config["virus_clustering"]["blast_threads"]
    shell:
        """
        # concatenate results
        cat {input.viruses} {input.virusdb} > {output.combined}

        # make a blast db from phage contigs
        makeblastdb -in {output.combined} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {output.combined} -db {params.blastdb} -out {params.blast_tsv} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90

        # calculate ani and af from blast results
        python {params.blastani_script} -i {params.blast_tsv} -o {params.ani_tsv}

        # cluster phage genomes based on 95% ani and 85% af
        python {params.cluster_script} --fna {output.combined} --ani {params.ani_tsv} --out {output.clusters} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


if config["virus_clustering"]["cluster_with_virusdb"]:
    clusters = (
        results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_clusters.tsv"
    )
else:
    clusters = results + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_clusters.tsv"


# extract cluster representatives from clusters
rule extract_cluster_representatives:
    input:
        viruses=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined.fna",
        clusters=clusters,
    output:
        results + "05_VIRUS_CLUSTERING/03_cluster_viruses/cluster_representatives.fna",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/05_extract_cluster_representatives.py.ipynb"


# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
rule virus_cluster_anlysis:
    input:
        clusters=clusters,
    output:
        results + "05_VIRUS_CLUSTERING/virus_clustering_figure.png",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/05_virus_clustering_analysis.py.ipynb"
