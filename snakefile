configfile: "config/config.yaml"


localrules:
    horesh_download,


rule all:
    input:
        expand("results/horesh/{name}", name=config["horesh_download"]),


rule horesh_download:
    output:
        "results/horesh/{name}",
    params:
        url=lambda w: config["horesh_download"][w.name],
    shell:
        "curl -L --fail --retry 5 --retry-delay 5 -o {output} '{params.url}'"
