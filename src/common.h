#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
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

template <typename T>
std::vector<std::size_t> sort_indices(const std::vector<T> &v) {
  std::vector<std::size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
       [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});
  return idx;
}

template <typename T>
double weighted_mean(const std::vector<T> &x, const std::vector<T> &w) {
    T weighted_sum = std::inner_product(x.begin(), x.end(), w.begin());
    T sum_of_weights = std::accumulate(w.begin(), w.end());
    return weighted_sum / sum_of_weights;
}

template <typename T>
double weighted_mean_product(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &w) {
    std::vector<T> product(x1.size());
    std::transform(x1.begin(), x1.end(), x2.begin(), x1.begin(), std::multiplies<T>());
    T weighted_sum = std::inner_product(product.begin(), product.end(), w.begin());
    T sum_of_weights = std::accumulate(w.begin(), w.end());
    return weighted_sum / sum_of_weights;
}

// template <typename T>
// double weighted_variance(const std::vector<T> &x, const std::vector<T> &w) {
//     T xmean = weighted_mean(x, w);
//     T weighted_sum = std::transform_reduce(
//         x.begin(),
//         x.end(),
//         w.begin(),
//         std::plus<T>(),
//         [=](const T& xval, const T& wval) {return wval * std::pow(xval - xmean, 2.0);}
//     );
//     T sum_of_weights = std::accumulate(w.begin(), w.end());
//     return weighted_sum / sum_of_weights;
// }

