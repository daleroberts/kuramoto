#pragma once
#include <Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

template <class Distribution>
class KuramotoProcess {
  public:
    KuramotoProcess(Vector x0, double dt, Distribution dist);
    
    Vector step();
    
  private:
    Vector _theta;
};