with open("assets/tr_options.txt", "r") as f:
    TR_LIST = [line.strip() for line in f]

rule all:
    input:
        "combined_outputs.txt"

rule compile_cpp:
    input: "src/cpp/peak_rp.cpp"
    output: "src/cpp/peak_rp"
    shell: "g++ -o {output} {input} -O3 -std=c++17 -pthread"

rule download_TR_target_sets:
    output:
        temp("data/target_sets/{distance}k/{tr}.tsv")
    params: 
        url="https://chip-atlas.dbcls.jp/data/hg19/target/{tr}.{distance}.tsv"
    log: 
        "logs/download_TR_target_sets/{distance}k/{tr}.log"
    shell:
        """
        wget -O {output} {params.url} 2> {log} || (echo "Failed to download {params.url}" >> {log}; touch {output})
        """

rule compile_TR_target_sets:
    input:
        lambda wildcards: expand(f"data/target_sets/{wildcards.distance}k/{{tr}}.tsv", tr=TR_LIST)
    output:
        "data/unweighted_network/{distance}k/ugrn.csv"
    log: 
        "logs/compile_TR_target_sets/{distance}k.log"
    shell:
        """
        python src/python/compile_TR_target_sets.py {input} {output} > {log} 2>&1
        """

rule calculate_H3K27ac_RP:
    input:
        "inputs/{context}.bed"
    output:
        "data/{context}/H3K27ac_RP.csv"
    params:
        refseq_genes="assets/hg19_refseq.tsv"
    shell:
        """
        src/cpp/peak_rp {input} {output} {params.refseq_genes}
        """

rule preprocess_TR_data:
    output:
        "preprocessed_TR_data.txt"
    shell:
        """
        # your code to download and preprocess TR data goes here
        """

rule intersect_TR_H3K27ac:
    input:
        "H3K27ac_data.txt",
        "preprocessed_TR_data.txt"
    output:
        "intersected_data.txt"
    shell:
        """
        # your code to intersect TR data with H3K27ac data goes here
        """

rule calculate_TR_RP:
    input:
        "intersected_data.txt"
    output:
        "TR_RP.txt"
    shell:
        """
        # your code to calculate TR RP goes here
        """

rule aggregate_TR_RPs:
    input:
        "TR_RP.txt"
    output:
        "aggregate_TR_RPs.txt"
    shell:
        """
        # your code to aggregate TR RPs goes here
        """

rule combine_outputs:
    input:
        "TR_target_sets.txt",
        "H3K27ac_RP.txt",
        "aggregate_TR_RPs.txt"
    output:
        "combined_outputs.txt"
    shell:
        """
        # your code to combine outputs from step 1, 2, and 4 goes here
        """
