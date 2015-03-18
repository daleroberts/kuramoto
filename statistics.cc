#include <algorithm>
#include <iomanip>
#include <limits>
#include <cmath>
#include "statistics.h"

Statistics::Statistics() {
    reset();
}

size_t Statistics::samples() const {
    return sample_number_;
}

double Statistics::weighted_sum() const {
    return sample_weight_;
}

double Statistics::mean() const {
    return sum_/sample_weight_;
}

double Statistics::variance() const {
    double m = mean();
    double v = quadratic_sum_/sample_weight_;
    v -= m*m;
    v *= sample_number_/(sample_number_-1.0);
    return v;
}

double Statistics::stddev() const {
    return std::sqrt(variance());
}

double Statistics::error_estimate() const {
    double var = variance();
    return std::sqrt(var/samples());
}

double Statistics::skewness() const {
    double s = stddev();

    if (s==0.0) return 0.0;

    double m = mean();
    double result = cubic_sum_/sample_weight_;
    result -= 3.0*m*(quadratic_sum_/sample_weight_);
    result += 2.0*m*m*m;
    result /= s*s*s;
    result *= sample_number_/(sample_number_-1.0);
    result *= sample_number_/(sample_number_-2.0);
    return result;
}

double Statistics::kurtosis() const {
    double m = mean();
    double v = variance();

    double c = (sample_number_-1.0)/(sample_number_-2.0);
    c *= (sample_number_-1.0)/(sample_number_-3.0);
    c *= 3.0;

    if (v==0) return c;

    double result = fourth_power_sum_/sample_weight_;
    result -= 4.0*m*(cubic_sum_/sample_weight_);
    result += 6.0*m*m*(quadratic_sum_/sample_weight_);
    result -= 3.0*m*m*m*m;
    result /= v*v;
    result *= sample_number_/(sample_number_-1.0);
    result *= sample_number_/(sample_number_-2.0);
    result *= (sample_number_+1.0)/(sample_number_-3.0);

    return result-c;
}

double Statistics::min() const {
    return min_;
}

double Statistics::max() const {
    return max_;
}

void Statistics::add(double value, double weight) {
    size_t old_samples = sample_number_;
    sample_number_++;

    sample_weight_ += weight;

    double temp = weight*value;
    sum_ += temp;
    temp *= value;
    quadratic_sum_ += temp;
    temp *= value;
    cubic_sum_ += temp;
    temp *= value;
    fourth_power_sum_ += temp;
    if (old_samples == 0) {
        min_ = max_ = value;
    } else {
        min_ = std::min(value, min_);
        max_ = std::max(value, max_);
    }
}

void Statistics::reset() {
    min_ = std::numeric_limits<double>::min();
    max_ = std::numeric_limits<double>::max();
    sample_number_ = 0;
    sample_weight_ = 0.0;
    sum_ = 0.0;
    quadratic_sum_ = 0.0;
    cubic_sum_ = 0.0;
    fourth_power_sum_ = 0.0;
}

Statistics Statistics::operator+(const Statistics& other) const {
    Statistics tmp = *this;
    tmp.sample_number_ += other.sample_number_;
    tmp.sample_weight_ += other.sample_weight_;
    tmp.sum_ += other.sum_;
    tmp.quadratic_sum_ += other.quadratic_sum_;
    tmp.cubic_sum_ += other.cubic_sum_;
    tmp.fourth_power_sum_ += other.fourth_power_sum_;
    tmp.min_ = std::min(other.min_, min_);
    tmp.max_ = std::max(other.max_, max_);
    return tmp;
}

Statistics& Statistics::operator=(const Statistics& other) {
    sample_number_ = other.sample_number_;
    sample_weight_ = other.sample_weight_;
    sum_ = other.sum_;
    quadratic_sum_ = other.quadratic_sum_;
    cubic_sum_ = other.cubic_sum_;
    fourth_power_sum_ = other.fourth_power_sum_;
    min_ = other.min_;
    max_ = other.max_;
    return *this;
}

Statistics& Statistics::operator+=(const Statistics& other) {
    *this = *this + other;
    return *this;
}
