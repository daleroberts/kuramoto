#pragma once
#include <cstdlib>

class Statistics {
  public:
    typedef double value_type;
    Statistics();

    size_t samples() const;
    double weighted_sum() const;
    double mean() const;
    double variance() const;
    double stddev() const;
    double error_estimate() const;
    double skewness() const;
    double kurtosis() const;
    double min() const;
    double max() const;
    void add(double value, double weight = 1.0);
    void reset();

    Statistics operator+(const Statistics& other);
    Statistics& operator=(const Statistics& other);
    Statistics& operator+=(const Statistics& other);

    template <class DataIterator> void add_from(DataIterator begin, DataIterator end) {
        for (;begin!=end;++begin)
            add(*begin);
    }

    template <class DataIterator, class WeightIterator>
    void add_from(DataIterator begin, DataIterator end, WeightIterator wbegin) {
        for (;begin!=end;++begin,++wbegin)
            add(*begin, *wbegin);
    }

  protected:
    size_t sample_number_;
    double sample_weight_;
    double sum_, quadratic_sum_, cubic_sum_, fourth_power_sum_;
    double min_, max_;
};
