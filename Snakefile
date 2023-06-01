# List of transcription regulators (TR) options
with open("assets/tr_options.txt", "r") as f:
    TR_LIST = [line.strip() for line in f]

# The rule 'all' is usually at the top and defines the final targets
rule all:
    input:
        "data/weighted_network/M.csv",
        "data/weighted_network/L.csv"

# Rule for compiling the cpp code
rule compile_cpp_code:
    input: "src/cpp/peak_rp.cpp"
    output: temp("src/cpp/peak_rp")
    shell: "g++ -o {output} {input} -O3 -std=c++17 -pthread"

# Rule for downloading transcription regulator target sets
rule download_target_sets:
    output: temp("data/target_sets/{tr}.tsv")
    params: 
        url="https://chip-atlas.dbcls.jp/data/hg19/target/{tr}.5.tsv"
    log: "logs/download_target_sets/{tr}.log"
    shell:
        """
        wget -O {output} {params.url} 2> {log} || (echo "Failed to download {params.url}" >> {log}; touch {output})
        """

# Rule for compiling transcription regulator target sets
rule compile_target_sets:
    input: expand("data/target_sets/{tr}.tsv", tr=TR_LIST)
    output: protected("data/unweighted_network/ugrn.csv")
    params: refseq_genes="assets/hg19_refseq.tsv"
    log: "logs/compile_target_sets.log"
    shell: "python src/python/compile_target_sets.py {input} {output} {params.refseq_genes} > {log} 2>&1"

# Rule for calculating H3K27ac RP
rule calculate_H3K27ac_RP:
    input:
        signal="inputs/{context}.bed",
        peak_rp="src/cpp/peak_rp"
    output: "data/{context}/H3K27ac_RP.csv"
    params: refseq_genes="assets/hg19_refseq.tsv"
    shell: "{input.peak_rp} {input.signal} {output} {params.refseq_genes}"

# Rule for downloading transcription regulator data
rule download_TR_data:
    output: temp("data/raw_TR_data/{tr}.bed")
    params: url="https://chip-atlas.dbcls.jp/data/hg19/assembled/Oth.ALL.50.{tr}.AllCell.bed"
    log: "logs/download_TR_data/{tr}.log"
    shell: "wget -O {output} {params.url} 2> {log} || echo 'Failed to download {params.url}' >> {log}"

# Rule for preprocessing transcription regulator data
rule preprocess_TR_data:
    input: "data/raw_TR_data/{tr}.bed"
    output: protected("data/TR_data/{tr}.bed")
    params: blacklist_file="assets/blacklists/ENCFF001TDO.bed"
    log: "logs/preprocess_TR_data/{tr}.log"
    shell: "bedtools intersect -v -wa -a {input} -b {params.blacklist_file} -sorted > {output}"

# Rule for intersecting transcription regulator data with H3K27ac
rule intersect_TR_H3K27ac:
    input: preprocessed_TR_data="data/TR_data/{tr}.bed"
    output: temp("data/{context}/intersected_data/{tr}.bed")
    params: h3K27ac_data="inputs/{context}.bed"
    log: "logs/intersect_TR_H3K27ac/{context}/{tr}.log"
    shell: "bedtools intersect -wa -a {input.preprocessed_TR_data} -b {params.h3K27ac_data} -sorted > {output}"

# Rule for binarizing chipseq binding
rule binarize_chipseq_binding:
    input: "data/{context}/intersected_data/{tr}.bed"
    output: temp("data/{context}/binarized_data/{tr}.bed")
    shell: "awk 'BEGIN {{FS=OFS=\"\t\"}} !seen[$1, $2, $3]++ {{print $1, $2, $3, 1}}' {input} > {output}"

# Rule for calculating transcription regulator RP
rule calculate_TR_RP:
    input:
        signal="data/{context}/binarized_data/{tr}.bed",
        peak_rp="src/cpp/peak_rp"
    output: temp("data/{context}/TR_RP/{tr}.csv")
    params: refseq_genes="assets/hg19_refseq.tsv"
    shell: "{input.peak_rp} {input.signal} {output} {params.refseq_genes}"

# Rule for aggregating transcription regulator RPs
rule aggregate_TR_RPs:
    input: lambda wildcards: expand(f"data/{wildcards.context}/TR_RP/{{tr}}.csv", tr=TR_LIST)
    output: "data/{context}/TR_RP.csv"
    shell: "python src/python/aggregate_TR_RPs.py {input} {output}"

# Rule for combining outputs into final CSV
rule combine_outputs:
    input:
        ugrn="data/unweighted_network/ugrn.csv",
        H3K27ac_RP="data/{context}/H3K27ac_RP.csv",
        TR_RP="data/{context}/TR_RP.csv"
    output: "data/weighted_network/{context}.csv"
    shell: "python src/python/combine_outputs.py {input.ugrn} {input.H3K27ac_RP} {input.TR_RP} {output}"
