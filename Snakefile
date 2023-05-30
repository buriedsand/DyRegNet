with open("assets/tr_options.txt", "r") as f:
    TR_LIST = [line.strip() for line in f]

rule all:
    input:
        UGRN5="data/unweighted_network/5k/ugrn.csv",
        M_H3K27ac_RP="data/M/H3K27ac_RP.csv",
        M_TR_RP="data/M/TR_RP.csv",
        L_H3K27ac_RP="data/L/H3K27ac_RP.csv",
        L_TR_RP="data/L/TR_RP.csv"


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
        protected("data/unweighted_network/{distance}k/ugrn.csv")
    log: 
        "logs/compile_TR_target_sets/{distance}k.log"
    shell:
        """
        python src/python/compile_TR_target_sets.py {input} {output} > {log} 2>&1
        """

rule calculate_H3K27ac_RP:
    input:
        signal="inputs/{context}.bed",
        peak_rp="src/cpp/peak_rp"
    output:
        "data/{context}/H3K27ac_RP.csv"
    params:
        refseq_genes="assets/hg19_refseq.tsv"
    shell:
        """
        {input.peak_rp} {input.signal} {output} {params.refseq_genes}
        """

rule download_TR_data:
    output:
        temp("data/raw_TR_data/{tr}.bed")
    params:
        url="https://chip-atlas.dbcls.jp/data/hg19/assembled/Oth.ALL.50.{tr}.AllCell.bed"
    log:
        "logs/download_TR_data/{tr}.log"
    shell:
        """
        wget -O {output} {params.url} 2> {log} || echo "Failed to download {params.url}" >> {log}
        """

rule preprocess_TR_data:
    input:
        "data/raw_TR_data/{tr}.bed"
    output:
        protected("data/TR_data/{tr}.bed")
    params:
        blacklist_file="assets/blacklists/ENCFF001TDO.bed"
    log: 
        "logs/preprocess_ChIP_data/{tr}.log"
    shell:
        """
        bedtools intersect -v -wa \
            -a {input} \
            -b {params.blacklist_file} \
            -sorted > {output}
        """

rule intersect_TR_H3K27ac:
    input:
        preprocessed_TR_data="data/TR_data/{tr}.bed"
    output:
        temp("data/{context}/intersected_data/{tr}.bed")
    params:
        h3K27ac_data="inputs/{context}.bed"
    log: 
        "logs/intersect_TR_H3K27ac/{context}/{tr}.log"
    shell:
        """
        bedtools intersect -wa \
            -a {input.preprocessed_TR_data} \
            -b {params.h3K27ac_data} \
            -sorted > {output}
        """

rule binarize_chipseq_binding:
    input:
        "data/{context}/intersected_data/{tr}.bed"
    output:
        temp("data/{context}/binarized_data/{tr}.bed")
    shell:
        """
        awk 'BEGIN {{FS=OFS=\"\t\"}} !seen[$1, $2, $3]++ {{print $1, $2, $3, 1}}' {input} > {output}
        """

rule calculate_TR_RP:
    input:
        signal="data/{context}/binarized_data/{tr}.bed",
        peak_rp="src/cpp/peak_rp"
    output:
        temp("data/{context}/TR_RP/{tr}.csv")
    params:
        refseq_genes="assets/hg19_refseq.tsv"
    shell:
        """
        {input.peak_rp} {input.signal} {output} {params.refseq_genes}
        """

rule aggregate_TR_RPs:
    input:
        lambda wildcards: expand(f"data/{wildcards.context}/TR_RP/{{tr}}.csv", tr=TR_LIST)
    output:
        "data/{context}/TR_RP.csv"
    shell:
        """
        python src/python/aggregate_TR_RPs.py {input} {output}
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
