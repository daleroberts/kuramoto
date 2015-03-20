#pragma once
#include <iostream>

class Statistics {
 public:
  Statistics();
  
 Statistics(size_t sample_number, double sample_weight, double sum,
	    double quadratic_sum, double cubic_sum, double fourth_power_sum,
	    double min, double max) : sample_number_(sample_number),
    sample_weight_(sample_weight), sum_(sum), quadratic_sum_(quadratic_sum),
    cubic_sum_(cubic_sum), fourth_power_sum_(fourth_power_sum), min_(min),
    max_(max) {}

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

  Statistics operator+(const Statistics& other) const; 
  Statistics& operator=(const Statistics& other);
  Statistics& operator+=(const Statistics& other);

  template <class DataIterator>
    void add_from(DataIterator begin, DataIterator end) {
    for (;begin!=end;++begin)
      add(*begin);
  }

  template <class DataIterator, class WeightIterator>
    void add_from(DataIterator begin, DataIterator end, WeightIterator wbegin) {
    for (;begin!=end;++begin,++wbegin)
      add(*begin, *wbegin);
  }

  template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & sample_number_;
      ar & sample_weight_;
      ar & sum_;
      ar & quadratic_sum_;
      ar & cubic_sum_;
      ar & fourth_power_sum_;
      ar & min_;
      ar & max_;
    }
    
  size_t sample_number_;
  double sample_weight_;
  double sum_, quadratic_sum_, cubic_sum_, fourth_power_sum_;
  double min_, max_;
};
