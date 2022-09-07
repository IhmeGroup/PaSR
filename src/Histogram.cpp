#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>

#include "common.h"
#include "Histogram.h"

Histogram::Histogram() :
    n_data(0), n_bins(0), s_rand(0)
{

}

void Histogram::clear() {
    n_bins = 0;
    n_data = 0;
    data.clear();
    bin_edges.clear();
}

void Histogram::printHist() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Histogram" << std::endl;
    std::cout << std::endl;
    std::cout <<
        std::left << std::setw(COL_WIDTH) << "Bin Start" <<
        std::left << std::setw(COL_WIDTH) << "Bin Stop" <<
        std::left << std::setw(COL_WIDTH) << "Count" << std::endl;
    for (int ib = 0; ib < n_bins; ib++) {
        std::cout <<
            std::left << std::setw(COL_WIDTH) << bin_edges[ib] << 
            std::left << std::setw(COL_WIDTH) << bin_edges[ib+1] <<
            std::left << std::setw(COL_WIDTH) << counts[ib] << std::endl;
    }
    std::cout << "--------------------------------------------------" << std::endl;
}

void Histogram::printPDF() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Histogram" << std::endl;
    std::cout << std::endl;
    std::cout <<
        std::left << std::setw(COL_WIDTH) << "Bin Start" <<
        std::left << std::setw(COL_WIDTH) << "Bin Stop" <<
        std::left << std::setw(COL_WIDTH) << "PDF" << std::endl;
    for (int ib = 0; ib < n_bins; ib++) {
        std::cout <<
            std::left << std::setw(COL_WIDTH) << bin_edges[ib] << 
            std::left << std::setw(COL_WIDTH) << bin_edges[ib+1] <<
            std::left << std::setw(COL_WIDTH) << pdf[ib] << std::endl;
    }
    std::cout << "--------------------------------------------------" << std::endl;

}

void Histogram::printCDF() {
    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "Histogram" << std::endl;
    std::cout << std::endl;
    std::cout <<
        std::left << std::setw(COL_WIDTH) << "Value" <<
        std::left << std::setw(COL_WIDTH) << "CDF" << std::endl;
    for (int ib = 0; ib < n_bins+1; ib++) {
        std::cout <<
            std::left << std::setw(COL_WIDTH) << bin_edges[ib] << 
            std::left << std::setw(COL_WIDTH) << cdf[ib] << std::endl;
    }
    std::cout << "--------------------------------------------------" << std::endl;

}

void Histogram::readHist(std::string hist_filename) {
    clear();
    std::string line, word;
    std::ifstream file(hist_filename);
    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream str(line);

            std::getline(str, word, ',');
            bin_edges.push_back(std::stod(word));

            if (std::getline(str, word, ',')) {
                counts.push_back(std::stod(word));
                n_bins++;
            }
        }
    } else {
        throw std::runtime_error("Histogram::readHist - could not open " + hist_filename + ".");
    }
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
    double bin_width = range() / (n_bins + 1);
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
    pdf.resize(n_bins);
    double total = 0.0;
// #pragma omp parallel for
    for (int ib = 0; ib < n_bins; ib++) {
        total += counts[ib] * (bin_edges[ib+1] - bin_edges[ib]);
    }
// #pragma omp parallel for
    for (int ib = 0; ib < n_bins; ib++) {
        pdf[ib] = counts[ib] / total;
    }
}

void Histogram::generateCDF() {
    cdf.resize(n_bins+1);
    cdf[0] = 0.0;
// #pragma omp parallel for
    for (int ib = 0; ib < n_bins; ib++) {
        cdf[ib+1] = cdf[ib] + (bin_edges[ib+1] - bin_edges[ib]) * pdf[ib];
    }
}

double Histogram::rand() {
    std::uniform_real_distribution<double> uni_real(0, 1);
    std::mt19937 eng(static_cast<uint64_t>(s_rand +
            std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())));
    s_rand++;
    return rand(uni_real, eng);
}

double Histogram::rand(std::uniform_real_distribution<double>& uni_real,
                       std::mt19937& eng) {
    double x_uni = uni_real(eng);
    for (int ib = 0; ib < n_bins; ib++) {
        if (x_uni < cdf[ib]) {
            continue;
        } else if (x_uni >= cdf[ib+1]) {
            continue;
        } else {
            return bin_edges[ib] + (x_uni - cdf[ib]) *
                (bin_edges[ib+1] - bin_edges[ib]) / (cdf[ib+1] - cdf[ib]);
        }
    }
    // Should never get here
    throw std::runtime_error("Histogram::rand - Invalid seed number. Distribution must have limits [0, 1).");
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

double Histogram::percentileToValue(double percentile) {
    if (percentile <= 0.0) {
        return bin_edges[0];
    } else if (percentile >= 1.0) {
        return bin_edges[n_bins+1];
    } else {
        for (int ib = 0; ib < n_bins; ib++) {
            if (percentile < cdf[ib]) {
                continue;
            } else if (percentile >= cdf[ib+1]) {
                continue;
            } else {
                return bin_edges[ib] + (percentile - cdf[ib]) *
                    (bin_edges[ib+1] - bin_edges[ib]) / (cdf[ib+1] - cdf[ib]);
            }
        }
    }
    // Should never get here
    throw std::runtime_error("Histogram::percentileToValue - Should never get here.");
}

double Histogram::valueToPercentile(double value) {
    if (value <= bin_edges[0]) {
        return 0.0;
    } else if (value >= bin_edges[n_bins+1]) {
        return 1.0;
    } else {
        for (int ib = 0; ib < n_bins; ib++) {
            if (value < bin_edges[ib]) {
                continue;
            } else if (value >= bin_edges[ib+1]) {
                continue;
            } else {
                return cdf[ib] + (value - bin_edges[ib]) *
                    (cdf[ib+1] - cdf[ib]) / (bin_edges[ib+1] - bin_edges[ib]);
            }
        }
    }
    // Should never get here
    throw std::runtime_error("Histogram::percentileToValue - Should never get here.");
}

Histogram::~Histogram() {

}
