# -----------------------------------------------------
# Virus identification (will always run)
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
# Virus identification rules
# -----------------------------------------------------
# -----------------------------------------------------
# 00 input
# -----------------------------------------------------
# symlink input paths to new paths
rule symlink_contigs:
    input:
        lambda wildcards: samples_df[
            (samples_df["group"] + samples_df["assembly"]) == wildcards.group_assembly
        ]["contigs"].iloc[0],
    output:
        results
        + "00_INPUT/{group_assembly}"
        + "_"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input} {output}
        """


# select which contig files to use depending on the files that are input
if len(config["virus_identification"]["assembled_contigs"]) == 0:
    contigs = (
        results
        + "03_READ_ASSEMBLY/04_contig_length_filter/{group_assembly}_"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    )
else:
    contigs = (
        results
        + "00_INPUT/{group_assembly}_"
        + config["read_assembly"]["assembly_output"]
        + ".fasta",
    )


# -----------------------------------------------------
# 01 MGV (& VirFinder)
# -----------------------------------------------------
# download build mgv repo and HMM files
rule build_mgv:
    output:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
    params:
        mgv_dir=resources + "mgv",
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm.gz",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm.gz",
    shell:
        """
        # clone MGV repository
        rm -rf {params.mgv_dir}
        git clone https://github.com/snayfach/MGV.git {params.mgv_dir}

        # download mgv hmm databases
        wget -O {params.imgvr_hmm} https://img.jgi.doe.gov//docs/final_list.hmms.gz
        wget -O {params.pfam_hmm} ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
        gunzip {params.imgvr_hmm}
        gunzip {params.pfam_hmm}
        """


# use prodigal to identify ORFs
rule mgv_prodigal:
    input:
        contigs,
    output:
        mgv_fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.fna",
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.faa",
        mgv_ffn=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.ffn",
    conda:
        "../envs/mgv.yml"
    shell:
        """        
        # call viral genes
        prodigal -i {input} \
        -a {output.mgv_faa} \
        -d {output.mgv_ffn}

        # create a symlink to original contigs fasta files
        ln -s {input} {output.mgv_fna}
        """


# use hmmsearch to find imgvr HMM hits
rule mgv_imgvr_hmmsearch:
    input:
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.faa",
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_imgvr.out",
    conda:
        "../envs/mgv.yml"
    threads: config["virus_identification"]["mgv_threads"]
    shell:
        """
        # identify imgvr viral genes using hmm search
        hmmsearch \
        -Z 1 \
        --cpu {threads} \
        --noali \
        --tblout {output} \
        {input.imgvr_hmm} {input.mgv_faa}
        """


# use hmmsearch to find pfam HMM hits
rule mgv_pfam_hmmsearch:
    input:
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.faa",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_pfam.out",
    conda:
        "../envs/mgv.yml"
    threads: config["virus_identification"]["mgv_threads"]
    shell:
        """
        # identify imgvr viral genes using hmmsearch
        hmmsearch \
        -Z 1 \
        --cpu {threads} \
        --noali \
        --tblout {output} \
        {input.pfam_hmm} {input.mgv_faa}
        """


# use mgv count_hmm_hits.py script to determine bacterial/viral HMM hits
rule mgv_count_hmm_hits:
    input:
        contigs=contigs,
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.faa",
        mgv_imgvr=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_imgvr.out",
        mgv_pfam=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_pfam.out",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_hmm_hits.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    conda:
        "../envs/mgv.yml"
    shell:
        """
        # change to mgv directory so scripts work correctly
        cd {params.mgv_dir}

        # count hmm hits using mgv script
        python count_hmm_hits.py \
        {input.contigs} {input.mgv_faa} \
        {input.mgv_imgvr} {input.mgv_pfam} > {output}
        """


# run virfinder to identify viral kmers
rule mgv_virfinder:
    input:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
        contigs=contigs,
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_virfinder.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    conda:
        "../envs/mgv.yml"
    shell:
        """
        # change to mgv directory so scripts work correclty
        cd {params.mgv_dir}

        # run mgv virfinder script
        Rscript virfinder.R \
        {input.contigs} {output}

        # convert NA to none
        sed -i 's/\tNA\tNA/\t0.0\t1.0/g' {output}
        """


# calculate the strand switch rate using mgv strand_switch.py script
rule mgv_strand_switch:
    input:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
        contigs=contigs,
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.faa",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_strand_switch.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    conda:
        "../envs/mgv.yml"
    shell:
        """
        # change to mgv directory so scripts work correctly
        cd {params.mgv_dir}

        # run mgv strand switch script
        python strand_switch.py \
        {input.contigs} {input.mgv_faa} > {output}
        """


# create master table using mgv master_table.py script
rule mgv_master_table:
    input:
        mgv_hmm_hits=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_hmm_hits.tsv",
        mgv_vf=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_virfinder.tsv",
        mgv_strand_switch=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_strand_switch.tsv",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_master_table.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    conda:
        "../envs/mgv.yml"
    shell:
        """
        # change to mgv directory so scripts run correctly
        cd {params.mgv_dir}

        # run mgv master table script
        python master_table.py \
        {input.mgv_hmm_hits} {input.mgv_vf} {input.mgv_strand_switch} > {output}
        """


# predict viral contigs using HMM hits, strand switch rate, and virfinder results (mgv viral_classify.py script)
rule mgv_viral_classify:
    input:
        mgv_fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}.fna",
        mgv_master_table=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_master_table.tsv",
    output:
        fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_final.fna",
        tsv=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_final.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
        input_base=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{group_assembly}",
        output_base=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_final",
    conda:
        "../envs/mgv.yml"
    shell:
        """
        # change to mgv directory so scripts run correctly
        cd {params.mgv_dir}

        # run mgv viral classify script
        python viral_classify.py \
        --features {input.mgv_master_table} \
        --in_base {params.input_base} \
        --out_base {params.output_base}
        """


# -----------------------------------------------------
# 02 VirSorter
# -----------------------------------------------------
# download virsorter db
rule download_virsorter_db:
    output:
        resources + "virsorter/virsorter-data-v2.tar.gz",
    params:
        output_dir=resources + "virsorter/",
    shell:
        """
        # download virsorter database
        wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz -P {params.output_dir}
        """


# extract virsorter db
rule extract_virsorter_db:
    input:
        resources + "virsorter/virsorter-data-v2.tar.gz",
    output:
        pfam=directory(resources + "virsorter/virsorter-data/PFAM_27/"),
        pgc=directory(resources + "virsorter/virsorter-data/Phage_gene_catalog/"),
        pgcpv=directory(
            resources + "virsorter/virsorter-data/Phage_gene_catalog_plus_viromes/"
        ),
        dockerfile=resources + "virsorter/virsorter-data/Dockerfile",
        ref_files=resources + "virsorter/virsorter-data/Generic_ref_file.refs",
        sup05=resources + "virsorter/virsorter-data/SUP05_SAGs_with_viruses.fna",
        vsrm=resources + "virsorter/virsorter-data/VirSorter_Readme.txt",
        vsrmv=resources + "virsorter/virsorter-data/VirSorter_Readme_viromes.txt",
    params:
        vs_dir=resources + "virsorter/",
    shell:
        """
        # download virsorter database
        tar -xvzf {input} -C {params.vs_dir}
        """


# identify viral contigs using virsorter
rule virsorter:
    input:
        contigs=contigs,
        vsrm=resources + "virsorter/virsorter-data/VirSorter_Readme.txt",
    output:
        vs_translation=results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{group_assembly}/fasta/input_sequences_id_translation.tsv",
        vs_results=results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{group_assembly}/Metric_files/VIRSorter_phage_signal.tab",
    params:
        output_dir=results + "04_VIRUS_IDENTIFICATION/02_virsorter/{group_assembly}",
        vs_dir=resources + "virsorter/virsorter-data/",
        extra_args=config["virus_identification"]["virsorter_arguments"],
    log:
        results + "00_LOGS/04_virus_identification_{group_assembly}.virsorter.log",
    conda:
        "../envs/virsorter.yml"
    threads: config["virus_identification"]["virsorter_threads"]
    shell:
        """
        rm -rf {params.output_dir}

        # run virsorter to identify viral contigs
        wrapper_phage_contigs_sorter_iPlant.pl \
        --fna {input.contigs} \
        --wdir {params.output_dir} \
        --ncpu {threads} \
        --data-dir {params.vs_dir} \
        {params.extra_args} > {log} 2>&1
        """


# -----------------------------------------------------
# 03 VirSorter2
# -----------------------------------------------------
# download virsorter2 db
rule download_virsorter2_db:
    output:
        resources + "virsorter2/Done_all_setup",
    params:
        vs2_dir=resources + "virsorter2/",
    conda:
        "../envs/virsorter2.yml"
    threads: config["virus_identification"]["virsorter2_threads"]
    shell:
        """
        # download virsorter2 database
        # remove the whole diretory specified by -d
        rm -rf db

        # run setup
        virsorter setup -d {params.vs2_dir} -j {threads}
        """


# run virsorter2 to identifiy viral contigs
rule virsorter2:
    input:
        contigs=contigs,
        vs2_db=resources + "virsorter2/Done_all_setup",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/03_virsorter2/{group_assembly}/final-viral-score.tsv",
    params:
        vs2_db=resources + "virsorter2",
        vs2_dir=results + "04_VIRUS_IDENTIFICATION/03_virsorter2/{group_assembly}",
        extra_args=config["virus_identification"]["virsorter2_arguments"],
    conda:
        "../envs/virsorter2.yml"
    threads: config["virus_identification"]["virsorter2_threads"]
    shell:
        """
        # run virsorter2
        virsorter run all --keep-original-seq \
        -w {params.vs2_dir} \
        -i {input.contigs} \
        -d {params.vs2_db} \
        -j {threads} \
        --prep-for-dramv \
        --min-score 0.0 \
        --rm-tmpdir \
        {params.extra_args}
        """


# -----------------------------------------------------
# 04 VIBRANT
# -----------------------------------------------------
# download vibrant database
rule download_vibrant_db:
    output:
        resources + "vibrant/db/databases/VIBRANT_setup.log",
    params:
        vb_dir=resources + "vibrant",
    conda:
        "../envs/vibrant.yml"
    shell:
        """
        rm -rf {params.vb_dir}

        # download vibrant database
        cd $CONDA_PREFIX/share/vibrant-1.2.1/db/databases
        ./VIBRANT_setup.py

        # create vibrant directory in resources folder
        mkdir {params.vb_dir}

        # move vibrant databases to resources folder
        cp -r $CONDA_PREFIX/share/vibrant-1.2.1/db {params.vb_dir}/db
        """


# run vibrant to identify viral contigs
rule vibrant:
    input:
        contigs=contigs,
        vb_db=resources + "vibrant/db/databases/VIBRANT_setup.log",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/04_vibrant/{group_assembly}/VIBRANT_{group_assembly}_contigs/VIBRANT_phages_{group_assembly}_contigs/{group_assembly}_contigs.phages_combined.txt",
    params:
        vb_db=resources + "vibrant/db/databases/",
        vb_files=resources + "vibrant/db/files/",
        vb_dir=results + "04_VIRUS_IDENTIFICATION/04_vibrant/{group_assembly}",
        extra_args=config["virus_identification"]["vibrant_arguments"],
    conda:
        "../envs/vibrant.yml"
    threads: config["virus_identification"]["vibrant_threads"]
    shell:
        """
        # remove vibrant directory
        rm -rf {params.vb_dir}

        # run vibrant for all virus types
        VIBRANT_run.py \
        -i {input.contigs} \
        -d {params.vb_db} \
        -m {params.vb_files} \
        -folder {params.vb_dir} \
        -t {threads} \
        {params.extra_args}
        """


# -----------------------------------------------------
# 05 DeepVirFinder
# -----------------------------------------------------
# build deepvirfinder to run
rule build_deepvirfinder:
    output:
        resources + "deepvirfinder/dvf.py",
    params:
        dvf_dir=resources + "/deepvirfinder",
    shell:
        """
        # git clone deepvirfinder repo
        git clone https://github.com/jessieren/DeepVirFinder {params.dvf_dir}
        """


# run deepvirfinder
rule deepvirfinder:
    input:
        contigs=contigs,
        dvf_script=resources + "deepvirfinder/dvf.py",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/05_deepvirfinder/{group_assembly}_contigs.fasta_gt1000bp_dvfpred.txt",
    params:
        output_dir=results + "04_VIRUS_IDENTIFICATION/05_deepvirfinder",
        model_dir=resources + "deepvirfinder/models",
    conda:
        "../envs/deepvirfinder.yml"
    threads: config["virus_identification"]["virfinder_threads"]
    shell:
        """
        # run deepvirfinder
        python {input.dvf_script} \
        -i {input.contigs} \
        -o {params.output_dir} \
        -l 1000 \
        -c {threads}
        """


# -----------------------------------------------------
# 06 Kraken2
# -----------------------------------------------------
# customize headers of virusdb for kraken
rule customize_virus_headers:
    input:
        genomes=config["virus_db"],
        metadata=config["virus_db_meta"],
    output:
        resources + "virusdb_kraken2db/kraken_formatted_viruses.fasta",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/customize_virus_headers.py.ipynb"


# build kraken database using custom virus database
rule kraken_build:
    input:
        resources + "virusdb_kraken2db/kraken_formatted_viruses.fasta",
    output:
        resources + "virusdb_kraken2db/hash.k2d",
    params:
        db=resources + "virusdb_kraken2db/",
    conda:
        "../envs/kraken2.yml"
    shell:
        """
        kraken2-build --download-taxonomy --db {params.db}
        kraken2-build --add-to-library {input} --db {params.db}
        kraken2-build --build --db {params.db}
        """


# align reads to kraken database
rule kraken2:
    input:
        db=resources + "virusdb_kraken2db/hash.k2d",
        contigs=contigs,
    output:
        classification=results
        + "04_VIRUS_IDENTIFICATION/06_kraken2/{group_assembly}/contigs.kraken2.classification.txt",
        report=results
        + "04_VIRUS_IDENTIFICATION/06_kraken2/{group_assembly}/kraken2.kreport",
    params:
        db=resources + "virusdb_kraken2db/",
    conda:
        "../envs/kraken2.yml"
    threads: config["virus_identification"]["kraken2_threads"]
    shell:
        """
        kraken2 {input.contigs} \
        --db {params.db} \
        --thread {threads} \
        --report {output.report} > {output.classification}
        """


# -----------------------------------------------------
# 07 Combine outputs
# -----------------------------------------------------
# determine input files for detecting virus sequences
if config["virus_identification"]["run_mgv"]:
    mgv1 = (
        results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_master_table.tsv",
    )
    mgv2 = (
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_final.tsv",
    )
else:
    mgv1 = contigs
    mgv2 = contigs
if config["virus_identification"]["run_virfinder"]:
    virfinder = (
        results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{group_assembly}_master_table.tsv",
    )
else:
    virfinder = contigs
if config["virus_identification"]["run_virsorter"]:
    virsorter = (
        results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{group_assembly}/Metric_files/VIRSorter_phage_signal.tab"
    )
    virsorter_translation = (
        results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{group_assembly}/fasta/input_sequences_id_translation.tsv"
    )
else:
    virsorter = contigs
    virsorter_translation = contigs
if config["virus_identification"]["run_virsorter2"]:
    virsorter2 = (
        results
        + "04_VIRUS_IDENTIFICATION/03_virsorter2/{group_assembly}/final-viral-score.tsv"
    )
else:
    virsorter2 = contigs
if config["virus_identification"]["run_vibrant"]:
    vibrant = (
        results
        + "04_VIRUS_IDENTIFICATION/04_vibrant/{group_assembly}/VIBRANT_{group_assembly}_contigs/VIBRANT_phages_{group_assembly}_contigs/{group_assembly}_contigs.phages_combined.txt",
    )
else:
    vibrant = contigs
if config["virus_identification"]["run_deepvirfinder"]:
    deepvirfinder = (
        results
        + "04_VIRUS_IDENTIFICATION/05_deepvirfinder/{group_assembly}_contigs.fasta_gt1000bp_dvfpred.txt",
    )
else:
    deepvirfinder = contigs
if config["virus_identification"]["run_kraken2"]:
    kraken2 = (
        results
        + "04_VIRUS_IDENTIFICATION/06_kraken2/{group_assembly}/contigs.kraken2.classification.txt",
    )
else:
    kraken2 = contigs


# combine outputs from all tool outputs
rule merge_reports_within_samples:
    input:
        contigs=contigs,
        mgv_results=mgv1,
        mgv_viruses=mgv2,
        vf_results=virfinder,
        vs_results=virsorter,
        vs_translation=virsorter_translation,
        vs2_results=virsorter2,
        vb_results=vibrant,
        dvf_results=deepvirfinder,
        kraken2_results=kraken2,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/07_combine_outputs/{group_assembly}/combined_report.csv",
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_virfinder=config["virus_identification"]["run_virfinder"],
        run_virsorter=config["virus_identification"]["run_virsorter"],
        run_virsorter2=config["virus_identification"]["run_virsorter2"],
        run_vibrant=config["virus_identification"]["run_vibrant"],
        run_deepvirfinder=config["virus_identification"]["run_deepvirfinder"],
        run_kraken2=config["virus_identification"]["run_kraken2"],
        group_assembly="{group_assembly}",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/04_merge_reports_within_samples.py.ipynb"


# combine viral contigs from all tool outputs using thresholds specified in config.yaml
rule merge_viral_contigs_within_samples:
    input:
        contigs=contigs,
        viral_report=results
        + "04_VIRUS_IDENTIFICATION/07_combine_outputs/{group_assembly}/combined_report.csv",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/07_combine_outputs/{group_assembly}/combined_viral_contigs.fasta",
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_vf=config["virus_identification"]["run_virfinder"],
        vf_score=config["virus_identification"]["virfinder_min_score"],
        run_vs=config["virus_identification"]["run_virsorter"],
        vs_cat=config["virus_identification"]["virsorter_cat"],
        run_vs2=config["virus_identification"]["run_virsorter2"],
        vs2_score=config["virus_identification"]["virsorter2_min_score"],
        run_dvf=config["virus_identification"]["run_deepvirfinder"],
        dvf_score=config["virus_identification"]["deepvirfinder_min_score"],
        run_vb=config["virus_identification"]["run_vibrant"],
        run_kraken2=config["virus_identification"]["run_kraken2"],
        assembly="{group_assembly}",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/04_merge_viral_contigs_within_samples.py.ipynb"


# -----------------------------------------------------
# Analyze combined virus data
# -----------------------------------------------------
# Determine inputs for combine_quast_outputs
def get_virus_identification_inputs(wildcards):
    assembly_report = pd.read_csv(
        checkpoints.combine_quast_outputs_across_samples.get(**wildcards).output[0],
        sep="\t",
    )
    return expand(
        results
        + "04_VIRUS_IDENTIFICATION/07_combine_outputs/{group_assembly}/combined_report.csv",
        group_assembly=list(set(assembly_report["Assembly"])),
    )


# combine virus reports
checkpoint combine_reports_across_samples:
    input:
        get_virus_identification_inputs,
    output:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# plot virus counts by tool
rule virus_identification_analysis:
    input:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    output:
        report(
            results + "04_VIRUS_IDENTIFICATION/virus_identification_figure.png",
            caption="../report/04_virus_identification_analysis.rst",
            category="Step 04: Virus identification",
        ),
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_vf=config["virus_identification"]["run_virfinder"],
        vf_score=config["virus_identification"]["virfinder_min_score"],
        run_vs=config["virus_identification"]["run_virsorter"],
        vs_cat=config["virus_identification"]["virsorter_cat"],
        run_vs2=config["virus_identification"]["run_virsorter2"],
        vs2_score=config["virus_identification"]["virsorter2_min_score"],
        run_dvf=config["virus_identification"]["run_deepvirfinder"],
        dvf_score=config["virus_identification"]["deepvirfinder_min_score"],
        run_vb=config["virus_identification"]["run_vibrant"],
        run_kraken2=config["virus_identification"]["run_kraken2"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/04_virus_identification_analysis.py.ipynb"
