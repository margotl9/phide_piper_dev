# -----------------------------------------------------
# Virus host (will always run)
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


# load sample information to be used in workflow
assemblies = list(set(samples_df["assembly"]))
contig_files = samples_df["contig"]


# -----------------------------------------------------
# Virus host rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 CRISPR Spacers
# -----------------------------------------------------
# make blastdb of uhgg spacers
rule make_spacer_blastdb:
    input:
        config["spacer_db"],
    output:
        resources + "crispr_spacers/crispr_spacers_db.ndb",
    params:
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db",
    conda:
        "../envs/blast.yml"
    shell:
        """
        # make blastdb
        makeblastdb \
        -in {input} \
        -out {params.spacers_blastdb} \
        -dbtype nucl
        """


# run crispropendb
rule blast_viruses_spacers:
    input:
        viruses=results
        + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db.ndb",
    output:
        results + "07_VIRUS_HOST/01_crispr_spacers/viruses_v_spacers_blast.tsv",
    params:
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db",
    conda:
        "../envs/blast.yml"
    threads: config["virus_host"]["blast_threads"]
    shell:
        """
        # blast viruses against uhgg spacers
        blastn \
        -query {input.viruses} \
        -db {params.spacers_blastdb} \
        -out {output} \
        -outfmt '6 std qlen slen' \
        -dust no \
        -word_size 18 \
        -num_threads {threads}
        """


# assign host taxonomy based on crispropendb blast
rule spacer_host_taxonomy:
    input:
        spacer_metadata=config["spacer_db_meta"],
        spacer_blast=results
        + "07_VIRUS_HOST/01_crispr_spacers/viruses_v_spacers_blast.tsv",
    output:
        report=results + "07_VIRUS_HOST/01_crispr_spacers/host_report.csv",
        taxonomy=results + "07_VIRUS_HOST/01_crispr_spacers/host_taxonomy.csv",
    params:
        max_mismatch=config["virus_host"]["max_spacer_mismatch"],
        min_spacer_coverage=config["virus_host"]["min_spacer_coverage"],
        min_agreement=config["virus_host"]["min_spacer_agreement"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/07_crispr_spacer_taxonomy.py.ipynb"


# -----------------------------------------------------
# 02 PHIST
# -----------------------------------------------------
# build phist
rule build_phist:
    output:
        resources + "phist/out/predictions.csv",
    params:
        phist_dir=resources + "phist",
    conda:
        "../envs/phist.yml"
    shell:
        """
        # git clone phist
        rm -rf {params.phist_dir}
        git clone --recurse-submodules https://github.com/refresh-bio/PHIST {params.phist_dir}

        # build phist
        cd {params.phist_dir}
        make
        mkdir ./out

        # test phist
        python3 phist.py ./example/virus ./example/host ./out/common_kmers.csv ./out/predictions.csv
        """


# organize viruses to by input into phist
rule split_viruses_for_phist:
    input:
        results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
    output:
        results + "07_VIRUS_HOST/02_phist/virus_fastas/viruses_prepared",
    params:
        fasta_dir=results + "07_VIRUS_HOST/02_phist/virus_fastas/",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/07_split_viruses_for_phist.py.ipynb"


# run phist using uhgg
rule phist:
    input:
        phist_build=resources + "phist/out/predictions.csv",
        virus_fastas=results + "07_VIRUS_HOST/02_phist/virus_fastas/viruses_prepared",
        bacteria=config["bacteria_db"],
    output:
        table=results + "07_VIRUS_HOST/02_phist/host_report.csv",
        predictions=results + "07_VIRUS_HOST/02_phist/host_predictions.csv",
    params:
        virus_dir=results + "07_VIRUS_HOST/02_phist/virus_fastas/",
        bacteria_db_dir=config["bacteria_db"].rpartition("/")[0],
        phist_script=resources + "phist/phist.py",
    threads: config["virus_host"]["phist_threads"]
    shell:
        """
        # run phist using uhgg
        python3 {params.phist_script} {params.virus_dir} {params.bacteria_db_dir} {output.table} {output.predictions} \
        --t {threads}
        """


# determine host taxonomy from refseq phist
rule phist_host_taxonomy:
    input:
        bacteria_db_metadata=config["bacteria_db_meta"],
        phist=results + "07_VIRUS_HOST/02_phist/host_predictions.csv",
    output:
        taxonomy=results + "07_VIRUS_HOST/02_phist/phist_host_taxonomy.csv",
        report=results + "07_VIRUS_HOST/02_phist/phist_host_report.csv",
    params:
        min_common_kmers=config["virus_host"]["min_phist_common_kmers"],
        min_agreement=config["virus_host"]["min_phist_agreement"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/07_phist_host_taxonomy.py.ipynb"


# -----------------------------------------------------
# 03 RAFAH
# -----------------------------------------------------
rule build_rafah:
    output:
        h3f=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f",
        h3i=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i",
        h3m=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m",
        h3p=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p",
        predict_script=resources + "rafah/RaFAH_Predict_Host.R",
        filename=resources + "rafah/MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData",
        rafah=resources + "rafah/RaFAH.pl",
        valid_cols=resources + "rafah/HP_Ranger_Model_3_Valid_Cols.txt",
    params:
        rafah_dir=resources + "rafah/",
    conda:
        "../envs/rafah.yml"
    shell:
        """
        # download rafah files
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/README.md
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_v0.3_hmm_models.tgz
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_v0.3_Ranger_Model.tgz
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH.pl
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/HP_Ranger_Model_3_Valid_Cols.txt
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_Predict_Host.R
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_Train_Model.R

        # decompress models
        cd {params.rafah_dir}
        tar -zxvf RaFAH_v0.3_hmm_models.tgz
        tar -zxvf RaFAH_v0.3_Ranger_Model.tgz
        """


rule rafah:
    input:
        h3f=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f",
        h3i=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i",
        h3m=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m",
        h3p=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p",
        rafah=resources + "rafah/RaFAH.pl",
        valid_cols=resources + "rafah/HP_Ranger_Model_3_Valid_Cols.txt",
        predict_script=resources + "rafah/RaFAH_Predict_Host.R",
        filename=resources + "rafah/MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData",
        viruses=results
        + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
    output:
        results + "07_VIRUS_HOST/03_rafah/Seq_Info_Prediction.tsv",
    params:
        virus_dir=results + "06_VIRUS_QUALITY/02_quality_filter/",
        hmm=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm",
        out_dir=results + "07_VIRUS_HOST/03_rafah",
    conda:
        "../envs/rafah.yml"
    threads: config["virus_host"]["rafah_threads"]
    shell:
        """
        cd {params.out_dir}

        perl {input.rafah} --predict \
        --genomes_dir {params.virus_dir} \
        --extension .fna \
        --valid_ogs_file {input.valid_cols} \
        --genomexog_table_file_name /home/carsonjm/CarsonJM/phide_piper/results/07_VIRUS_HOST/03_rafah/RaFAH_Genome_to_OGs_Score_Min_Score_50-Max_evalue_1e-05_Prediction.tsv \
        --hmmer_db_file_name {params.hmm} \
        --r_script_predict_file_name {input.predict_script} \
        --r_model_file_name {input.filename} \
        --threads {threads}
        """


# -----------------------------------------------------
# 04 HostG
# -----------------------------------------------------
# rule build_hostg:
#     output:

#     params:
#         hostg_dir=resources + "hostg",
#         dataset=resources + "hostg/dataset",
#     conda:
#         "../envs/hostg.yml"
#     threads: config["virus_host"]["hostg_threads"]
#     shell:
#         """

#         git clone https://github.com/KennthShang/HostG.git {params.hostg_dir}

#         cd {params.dataset}
#         bzip2 -d protein.fasta.bz2
#         bzip2 -d nucl.fasta.bz2

#         """


# rule hostg:
#     input:
#         script=resources + "hostg/run_Speed_up.py",
#         viruses=results
#         + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
#     output:
#         results + "07_VIRUS_HOST/04_rafah/"
#     params:
#         out_dir=results + "07_VIRUS_HOST/04_rafah/final_prediction.csv"
#     conda:
#     threads:
#     shell:
#         """
#         cd {params.out_dir}
#         python {input.script} \
#         --contigs {input.viruses} \
#         --len 1000 \
#         --t 0
#         """
