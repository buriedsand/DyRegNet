#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <sstream>

// Import necessary libraries for handling CSV files
#include "csv.h" // Use a CSV parsing library, e.g., https://github.com/ben-strasser/fast-cpp-csv-parser

const int L = 3000; // 100000;
const int delta = 300; // 10000;
const double mu = -std::log(static_cast<double>(L) / (3 * delta));

double calculate_single_gene_rp(const std::vector<std::tuple<int, int, double>> &sample_values, int t_k) {
    double R_jk = 0;
    for (const auto &[start, end, s_ij] : sample_values) {
        int adjusted_start = std::max(start, t_k - L);
        int adjusted_end = std::min(end, t_k + L);
        for (int i = adjusted_start; i < adjusted_end; ++i) {
            double d = std::abs(i - t_k) / static_cast<double>(L);
            double exp_neg_mu_d = std::exp(-mu * d);
            double w_i = 2 * exp_neg_mu_d / (1 + exp_neg_mu_d);
            R_jk += w_i * s_ij;
        }
    }
    return R_jk;
}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <refseq_file>" << std::endl;
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];
    const char *refseq_file = argv[3];

    // Read signal values (e.g., binarized TF binding, H3K27ac signal)
    io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>> signal_values_reader(input_file);

    std::string chr;
    int start, end;
    double s_ij;
    std::unordered_map<std::string, std::vector<std::tuple<int, int, double>>> signal_values;

    while (signal_values_reader.read_row(chr, start, end, s_ij)) {
        signal_values[chr].push_back({start, end, s_ij});
    }

    // Read refseq genes
    io::CSVReader<4, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>> refseq_genes_reader(refseq_file);

    std::string gene_id, strand;
    int tss;
    std::vector<std::tuple<std::string, std::string, int>> refseq_genes;

    while (refseq_genes_reader.read_row(gene_id, chr, strand, tss)) {
        refseq_genes.push_back({gene_id, chr, tss});
    }


    // Calculate peak RP
    std::unordered_map<std::string, double> peak_rp;

    for (const auto &[gene_id, gene_chr, t_k] : refseq_genes) {
        int interval_start = std::max(t_k - L, 0);
        int interval_end = t_k + L;

        std::vector<std::tuple<int, int, double>> interval_values;
        for (const auto &row : signal_values[gene_chr]) {
            int start, end;
            double s_ij;
            std::tie(start, end, s_ij) = row;
            if (end >= interval_start && start <= interval_end) {
                interval_values.push_back(row);
            }
        }


        double R_jk = calculate_single_gene_rp(interval_values, t_k);
        peak_rp[gene_id] = R_jk;
    }


    // Write the output to a file
    std::ofstream outfile(output_file);

    if (!outfile.is_open()) {
        std::cerr << "Error opening output file: " << output_file << std::endl;
        return 1;
    }

    for (const auto &[key, value] : peak_rp) {
        outfile << key << "," << value << std::endl;
    }

    outfile.close();
    return 0;
}

   
