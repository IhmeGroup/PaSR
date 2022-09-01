#include <iostream>

#include "cantera/base/ctexceptions.h"

class Histogram {
public:
    explicit Histogram();
    ~Histogram();

    void clear();
    void printHist();
    void printPDF();
    void printCDF();

    void readData(std::string data_filename);
    void resizeData(int n_data_);
    void setData(const double* data_);

    void readHist(std::string hist_filename);
    void generate(int n_bins_=0);
    void generateHist(int n_bins_=0);
    void generatePDF();
    void generateCDF();

    double rand();
    double rand(std::uniform_real_distribution<double>& uni_real,
                std::mt19937& eng);
    double mean();
    double median();
    double mode();
    double min();
    double max();
    double range();
    double percentileToValue(double percentile);
    double valueToPercentile(double value);
    double variance();
    double stddev();
    double skewness();
    double kurtosis();

    double* getBinEdges() {
        return bin_edges.data();
    }

    double* getCounts() {
        return counts.data();
    }

    int dataSize() {
        return n_data;
    }

    int nBins() {
        return n_bins;
    }

    double* getPDF() {
        return pdf.data();
    }
    
    double* getCDF() {
        return cdf.data();
    }

protected:
    int n_data;
    int n_bins;
    std::vector<double> data;
    std::vector<double> bin_edges;
    std::vector<double> counts;
    std::vector<double> pdf;
    std::vector<double> cdf;

    int s_rand = 0;

private:

};