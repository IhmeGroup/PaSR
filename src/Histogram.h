#include <iostream>

#include "cantera/base/ctexceptions.h"

class Histogram {
public:
    explicit Histogram();
    ~Histogram();

    void read(std::string hist_filename);
    void resizeData(int n_data_);
    void setData(const double* data_);

    void generateHist(int n_bins_=0);
    void generatePDF();
    void generateCDF();

    double mean();
    double median();
    double mode();
    double min();
    double max();
    double range();
    double variance();
    double stddev();
    double skewness();
    double kurtosis();

    double rand();
    double* pdf();
    double* cdf();

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

protected:
    int n_data;
    int n_bins;
    std::vector<double> data;
    std::vector<double> bin_edges;
    std::vector<double> counts;
    std::vector<double> pdf;
    std::vector<double> cdf;

private:

};