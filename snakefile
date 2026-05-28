configfile: "config/config.yaml"


localrules:
    horesh_download,
    extract_focal_cons_seq,
    extract_focal_gene_info,
    parse_homology,


rule all:
    input:
        expand("results/horesh/{name}", name=config["horesh_download"]),
        "results/homology/focal_vs_horesh_mapping.csv",


rule horesh_download:
    output:
        "results/horesh/{name}",
    params:
        url=lambda w: config["horesh_download"][w.name],
    shell:
        "curl -L --fail --retry 5 --retry-delay 5 -o {output} '{params.url}'"


rule extract_focal_cons_seq:
    input:
        config["focal_genes"]["cons_seq"],
    output:
        "results/focal_genes/consensus_sequences.fa",
    shell:
        "mkdir -p $(dirname {output}) && tar -xJf {input} -C $(dirname {output})"


rule extract_focal_gene_info:
    input:
        config["focal_genes"]["gene_info"],
    output:
        "results/focal_genes/gene_info.csv",
    shell:
        "mkdir -p $(dirname {output}) && tar -xJf {input} -C $(dirname {output})"


rule blastn_focal_genes:
    input:
        query=rules.extract_focal_cons_seq.output,
        db_fa="results/horesh/F3_pan_genome_reference.fa",
    output:
        tsv="results/blastn/focal_vs_horesh_raw.tsv",
    log:
        "logs/blastn/focal_vs_horesh.log",
    conda:
        "config/conda_envs/blast.yaml"
    threads: 4
    params:
        db="results/blastn/db",
        min_pident=config["homology"]["min_pident"],
    shell:
        """
        exec >{log} 2>&1
        mkdir -p $(dirname {params.db})
        makeblastdb -in {input.db_fa} -dbtype nucl -out {params.db}

        blastn \
            -query {input.query} \
            -db {params.db} \
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
            -num_threads {threads} \
            -evalue 1e-5 \
            -perc_identity {params.min_pident} \
            -max_target_seqs 20 \
            -out {output.tsv}
        """


rule parse_homology:
    input:
        hits=rules.blastn_focal_genes.output.tsv,
        query_fa=rules.extract_focal_cons_seq.output,
    output:
        mapping="results/homology/focal_vs_horesh_mapping.csv",
        unmatched="results/homology/focal_vs_horesh_unmatched.txt",
    log:
        "logs/homology/focal_vs_horesh.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    params:
        min_pident=config["homology"]["min_pident"],
        min_qcovs=config["homology"]["min_qcovs"],
        min_scovs=config["homology"]["min_scovs"],
    shell:
        """
        python3 scripts/parse_homology.py \
            --hits {input.hits} \
            --query-fa {input.query_fa} \
            --mapping {output.mapping} \
            --unmatched {output.unmatched} \
            --min-pident {params.min_pident} \
            --min-qcovs {params.min_qcovs} \
            --min-scovs {params.min_scovs} \
            >{log} 2>&1
        """
