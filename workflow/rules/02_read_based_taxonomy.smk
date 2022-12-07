# -------------------------------------
# Read-Based Taxonomy Module
# -------------------------------------
import pandas as pd
import os

# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv("config/samples.tsv", sep="\t")

# get current working directory so absolute paths can be used for input/output files
results = os.getcwd()

# load report
report: "../report/workflow.rst"

# load resources folder path
resources = config["resources_path"]

# load sample information to be used in workflow
samples = samples_df["sample"]
samples_assemblies=list(set(samples_df["sample"].astype("string") + "_" + samples_df["assembly"].astype("string")))


# -------------------------------------
# Read-Based Taxonomy Rules
# -------------------------------------
### MetaPhlan ###
# download metaphlan database
rule metaphlan_database:
    output:
        metphlan_db=resources + "/metaphlan/mpa_v30_CHOCOPhlAn_201901/mpa_v30_CHOCOPhlAn_201901.1.bt2",
    params:
        metaphlan_db=resources + "/metaphlan/mpa_v30_CHOCOPhlAn_201901/",
    conda: 
        "../envs/metaphlan.yml"
    shell:
        """
        # install metaphlan database
        metaphlan --install \
        --index mpa_v30_CHOCOPhlAn_201901 \
        --bowtie2db {params.metaphlan_db}
        """


# taxonomically annotate reads with metaphlan
rule metaphlan:
    input:
        metphlan_db=resources + "/metaphlan/mpa_v30_CHOCOPhlAn_201901/mpa_v30_CHOCOPhlAn_201901.1.bt2",
        R1=results + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_1.fastq.gz",
        R2=results + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_paired_2.fastq.gz",
        R1S=results + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_1.fastq.gz",
        R2S=results + "/01_READ_PREPROCESSING/03_kneaddata/{sample_assembly}_unmatched_2.fastq.gz",
    output:
        profile=results+"/02_READ_BASED_TAXONOMY/01_metaphlan/{sample_assembly}_profile.txt",
        bowtie2=results+"/02_READ_BASED_TAXONOMY/01_metaphlan/{sample_assembly}_bowtie2.bz2",
    params:
        merged_fastq=results+"/02_READ_BASED_TAXONOMY/01_metaphlan/{sample_assembly}_merged.fastq",
        metaphlan_db=resources + "/metaphlan/mpa_v30_CHOCOPhlAn_201901/",
        extra_args=config['metaphlan']['extra_args'],
    threads:
        10
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        # merge all fastq files from a sample
        zcat {input.R1} {input.R2} {input.R1S} {input.R2S} > {params.merged_fastq}
        
        # run metaphlan to taxonomically annotate reads
        metaphlan {params.merged_fastq} --output_file {output.profile} \
        --bowtie2db {params.metaphlan_db} \
        --bowtie2out {output.bowtie2} \
        --nproc {threads} \
        {params.extra_args}

        # remove merged fastq files
        rm {params.merged_fastq}
        """


# combine metaphlan results
rule combine_metaphlan_profiles:
    input:
        expand(results+"/02_READ_BASED_TAXONOMY/01_metaphlan/{sample_assembly}_profile.txt", sample_assembly=samples_assemblies)
    output:
        results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles.txt",
    params:
        profiles_dir=results+"/02_READ_BASED_TAXONOMY/01_metaphlan/",
    conda:
        "../envs/metaphlan.yml"
    shell:
        """
        # merge metaphlan results
        merge_metaphlan_tables.py \
        {input} > {output}
        """

# visualize metaphlan results
rule metaphlan_level_reduction:
    input:
        results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles.txt"
    output:
        species=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_species.csv",
        genus=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_genus.csv",
        family=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_family.csv",
        order=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_order.csv",
        class_df=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_class.csv",
        phylum=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_phylum.csv",
        kingdom=results+"/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_kingdom.csv",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/02_read_based_taxonomy_metaphlan_levels.py"


# visualize metaphlan results
rule metaphlan_analysis:
    input:
        results + "/02_READ_BASED_TAXONOMY/metaphlan_combined_profiles_phylum.csv",
    output:
        figure=report(
            results + "/02_READ_BASED_TAXONOMY/metaphlan_phyla_barplots.png",
            caption="../report/02_read_based_taxonomy_metaphlan.rst",
            category="Step 02: Read-based taxonomy",
        ),
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/02_read_based_taxonomy_metaphlan_analysis.py"