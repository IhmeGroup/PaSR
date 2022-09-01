#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <omp.h>

#define COL_WIDTH 15

#pragma omp declare \
    reduction( \
        vec_double_min : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            [](double x1, double x2) { return std::min(x1, x2); })) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare \
    reduction( \
        vec_double_max : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            [](double x1, double x2) { return std::max(x1, x2); })) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))
#pragma omp declare \
    reduction( \
        vec_double_plus : \
        std::vector<double> : \
        std::transform( \
            omp_out.begin(), \
            omp_out.end(), \
            omp_in.begin(), \
            omp_out.begin(), \
            std::plus<double>())) \
    initializer( \
        omp_priv = decltype(omp_orig)(omp_orig.size()))