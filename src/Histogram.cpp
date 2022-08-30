#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>

#include "common.h"
#include "Histogram.h"

Histogram::Histogram() {

}

void Histogram::generate(int n_bins_) {
    generateHist(n_bins_);
    generatePDF();
    generateCDF();
}

void Histogram::generateHist(int n_bins_) {

    // Determine bin count if necessary
    if (n_bins_ > 0) {
        n_bins = n_bins_;
    } else {
        // Sturges' formula
        n_bins = 1 + std::ceil(std::log2(n_data));
    }

    bin_edges.resize(n_bins+1);
    counts.resize(n_bins);

    // Min and max element are at midpoint of first and last bins
    bin_width = range() / (n_bins + 1);
    double bin_start = min() - 0.5*bin_width;

    // Fill bin edges
    for (int ib = 0; ib < n_bins+1; ib++) {
        bin_edges[ib] = bin_start + ib * bin_width;
    }

    // Compute counts
    std::vector<double> counts_temp(n_bins, 0.0);
#pragma omp parallel for reduction(vec_double_plus:counts_temp)
    for (int id = 0; id < n_data; id++) {
        for (int ib = 0; ib < n_bins; ib++) {
            if (data[id] < bin_edges[ib]) {
                continue;
            } else if (data[id] >= bin_edges[ib+1]) {
                continue;
            } else {
                counts_temp[ib]++;
            }
        }
    }

    // Fill counts
    for (int ib = 0; ib < n_bins; ib++) {
        counts[ib] = counts_temp[ib];
    }
}

void Histogram::generatePDF() {
    double total = bin_width * n_data;
// #pragma omp parallel for
    for (int ib = 0; ib < n_bins; ib++) {
        pdf[ib] = counts[ib] / total;
    }
}

void Histogram::generateCDF() {
    cdf[0] = 0.0;
    for (int ib = 0; ib < n_bins-1; ib++) {
        cdf[ib+1] = cdf[ib] + bin_width * pdf[ib];
    }
}

double Histogram::rand(std::uniform_real_distribution<double>& uni_real,
                       std::mt19937& eng) {
    double x_uni = uni_real(eng);
    return x_uni;
}

double Histogram::min() {
    return *std::min_element(data.begin(), data.end());
}

double Histogram::max() {
    return *std::max_element(data.begin(), data.end());
}

double Histogram::range() {
    return max() - min();
}

Histogram::~Histogram() {

}