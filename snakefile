configfile: "config/config.yaml"


localrules:
    horesh_download,


rule all:
    input:
        expand("results/horesh/{name}", name=config["horesh_download"]),
        "results/presence_absence/focal_vs_horesh.csv",
        "results/figs/st_boxplots.png",
        "results/figs/st_boxplots_with_counts.png",
        "results/figs/gene_enrichment_scatter.png",
        "results/figs/st_isolates_per_st.png",
        "results/figs/st_pairwise_distance.png",
        "results/figs/phylo_heatmap.png",


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
        "logs/blastn_focal_genes.log",
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
        "logs/parse_homology.log",
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


rule focal_genes_presence_absence:
    input:
        mapping=rules.parse_homology.output.mapping,
        pa="results/horesh/F4_complete_presence_absence.csv",
    output:
        "results/presence_absence/focal_vs_horesh.csv",
    log:
        "logs/focal_genes_presence_absence.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    shell:
        """
        python3 scripts/focal_genes_presence_absence.py \
            --mapping {input.mapping} \
            --presence-absence {input.pa} \
            --output {output} \
            >{log} 2>&1
        """


rule plot_st_boxplots:
    input:
        pa=rules.focal_genes_presence_absence.output,
        gene_info=rules.extract_focal_gene_info.output,
        metadata="results/horesh/F1_genome_metadata.csv",
    output:
        boxplot=multiext("results/figs/st_boxplots", ".png", ".pdf"),
        with_counts=multiext("results/figs/st_boxplots_with_counts", ".png", ".pdf"),
    log:
        "logs/plot_st_boxplots.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    params:
        boxplot_base="results/figs/st_boxplots",
        with_counts_base="results/figs/st_boxplots_with_counts",
        high_load_threshold=config["high_load_threshold"],
    shell:
        """
        python3 scripts/plots/plot_st_boxplots.py \
            --pa {input.pa} \
            --gene-info {input.gene_info} \
            --metadata {input.metadata} \
            --out-boxplot-base {params.boxplot_base} \
            --out-boxplot-with-counts-base {params.with_counts_base} \
            --high-load-threshold {params.high_load_threshold} \
            >{log} 2>&1
        """


rule plot_gene_enrichment_scatter:
    input:
        pa=rules.focal_genes_presence_absence.output,
        gene_info=rules.extract_focal_gene_info.output,
        metadata="results/horesh/F1_genome_metadata.csv",
        highlight="config/focal_genes/genes_to_highlight.csv",
    output:
        fig=multiext("results/figs/gene_enrichment_scatter", ".png", ".pdf"),
    log:
        "logs/plot_gene_enrichment_scatter.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    params:
        out_base="results/figs/gene_enrichment_scatter",
        high_load_threshold=config["high_load_threshold"],
    shell:
        """
        python3 scripts/plots/plot_gene_enrichment_scatter.py \
            --pa {input.pa} \
            --gene-info {input.gene_info} \
            --metadata {input.metadata} \
            --highlight {input.highlight} \
            --out-fig-base {params.out_base} \
            --high-load-threshold {params.high_load_threshold} \
            >{log} 2>&1
        """


rule plot_st_phylo_distance:
    input:
        tree="results/horesh/tree_500.nwk",
        metadata="results/horesh/F1_genome_metadata.csv",
    output:
        isolates=multiext("results/figs/st_isolates_per_st", ".png", ".pdf"),
        distance=multiext("results/figs/st_pairwise_distance", ".png", ".pdf"),
        counts_csv="results/figs/st_isolate_counts.csv",
        distance_csv="results/figs/st_pairwise_distance.csv",
    log:
        "logs/plot_st_phylo_distance.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    params:
        isolates_base="results/figs/st_isolates_per_st",
        distance_base="results/figs/st_pairwise_distance",
        tree_base="results/figs/tree_phylogroup",
    shell:
        """
        python3 scripts/plots/plot_st_phylo_distance.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --out-isolates-base {params.isolates_base} \
            --out-distance-base {params.distance_base} \
            --out-tree-base {params.tree_base} \
            --out-counts-csv {output.counts_csv} \
            --out-distance-csv {output.distance_csv} \
            >{log} 2>&1
        """


rule plot_phylo_heatmap:
    input:
        pa=rules.focal_genes_presence_absence.output,
        gene_info=rules.extract_focal_gene_info.output,
        metadata="results/horesh/F1_genome_metadata.csv",
        tree="results/horesh/tree_500.nwk",
    output:
        fig=multiext("results/figs/phylo_heatmap", ".png", ".pdf"),
        pa_csv="results/figs/phylo_heatmap_pa_matrix.csv",
    log:
        "logs/plot_phylo_heatmap.log",
    conda:
        "config/conda_envs/bioinfo.yml"
    params:
        out_base="results/figs/phylo_heatmap",
    shell:
        """
        python3 scripts/plots/plot_phylo_heatmap.py \
            --pa {input.pa} \
            --gene-info {input.gene_info} \
            --metadata {input.metadata} \
            --tree {input.tree} \
            --out-heatmap-base {params.out_base} \
            --out-pa-csv {output.pa_csv} \
            >{log} 2>&1
        """
