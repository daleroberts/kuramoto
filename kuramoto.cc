/* Kuramoto Model with Tempered Stable Noise
 * Dale Roberts <dale.o.roberts@gmail.com>
 */

#include <random>
#include <cmath>
#include <Eigen/Dense>
#include <stdio.h>
#include <iostream>
#include <unordered_map>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

class Graph {
public:
    Graph(const std::string& str) {
        std::string s(str);
        if (s.find(">>graph6<<") != std::string::npos)
            s = s.substr(10);

        auto data = graph6_to_data(s);

        size_t n;
        if (data[0] <= 62) {
            n = data[0];
            data.erase(data.begin());
        } else if (data[1] <= 62) {
            n = (data[1]<<12) + (data[2]<<6) + data[3];
            data.erase(data.begin(), data.begin() + 4);
        } else {
            n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18)
                + (data[5]<<12) + (data[6]<<6) + data[7];
            data.erase(data.begin(), data.begin() + 8);
        }
        n_ = n;

        size_t nd = (n*(n-1)/2 + 5) / 6;
        if (data.size() != nd)
            std::cout << "error parsing graph: " << s << std::endl;

        std::vector<size_t> nodes;
        for (size_t i = 0; i < n; ++i)
            nodes.push_back(i);

        size_t k = 5;
        auto d = data.begin();
        for (size_t j = 1; j < n; ++j) {
            for (size_t i = 0; i < j; ++i) {
                if ((*d>>k)&1)
                    add_edge(i, j);
                if (!k--) {
                    k = 5;
                    d++;
                }
            }
        }
    }

    size_t size() {
        return n_;
    }

    void add_nodes_from(std::vector<size_t> nodes) {
        for (auto &n : nodes) {
            if (!has_node(n)) {
                nodemap m;
                adj_[n] = m;
            }
        }
    }

    std::vector<size_t> nodes(void) {
        std::vector<size_t> v;
        for (auto &el : nodes_)
            v.push_back(el.first);
        std::sort(v.begin(), v.end());
        return v;
    }

    void add_edge(const size_t i, const size_t j) {
        if (!has_node(i)) {
            adj_[i];
            nodes_[i];
        }
        if (!has_node(j)) {
            adj_[j];
            nodes_[j];
        }
        adj_[i][j] = true; 
        adj_[j][i] = true; 
    }

    bool has_node(const size_t i) {
        return nodes_.count(i) > 0;
    }

    std::vector<size_t> neighbours(const size_t i) {
        std::vector<size_t> v;
        for (auto &a : adj_[i])
            v.push_back(a.first);
        std::sort(v.begin(), v.end());
        return v;
    }

private:
    typedef std::unordered_map<size_t, bool> nodemap;
    nodemap nodes_;
    std::unordered_map<size_t, nodemap> adj_;
    size_t n_;

    std::vector<int> graph6_to_data(std::string s) {
        std::vector<int> v;
        for (auto &c : s)
            v.push_back((int)c-63);
        int max = *std::max_element(v.begin(), v.end());
        int min = *std::min_element(v.begin(), v.end());
        if (v.size() > 0 && (min < 0 || max > 63))
            v.clear();
        return v;
    }
};

using namespace std;

inline double mod2pi(double theta) {
    while (theta > 2*M_PI)
        theta -= 2*M_PI;
    while (theta < 0.0)
        theta += 2*M_PI;
    return theta;
}

template <typename Derived>
inline double order_param(Eigen::MatrixBase<Derived>& theta) {
    int N = theta.size();

    double sum_real = 0.0;
    double sum_complex = 0.0;
    for (int j = 0; j < N; j++) {
        sum_real += cos(theta(j));
        sum_complex += sin(theta(j));
    }

    return sqrt(sum_real*sum_real + sum_complex*sum_complex)/N;
}


inline double sinc(double x) {
    double ax = fabs(x);
    if (ax < 0.006) {
        if (x == 0.)
            return 1;
        double x2 = x*x;
        if(ax < 2e-4)
            return 1. - x2/6.;
        else 
            return 1. - x2/6.*(1 - x2/20.);
    }
    /* else */
    return sin(x)/x;
}

inline double A(double x, double rho) {
    double Irho = 1.-rho;
    return pow(Irho*sinc(Irho*x),Irho)*pow(rho*sinc(rho*x),rho)/sinc(x);
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
void path(const double globalCoupling,
        const double noiseScale,
        Eigen::MatrixBase<DerivedA>& initialConditions,
        Eigen::MatrixBase<DerivedB>& intrinsicFreqs,
        Eigen::MatrixBase<DerivedC>& adjacencyMatrix,
        const double rho,
        const double c,
        const double alpha,
        const double p,
        const double t,
        const int timeSteps,
        const int seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> runif(0,1);
    std::exponential_distribution<> rexp(1);

    double dt = t/timeSteps;

    double gamma = c/(1-rho) - p;
    double cf = pow(cos(M_PI_2*rho),-1./rho);
    double sigma = pow(-dt*c*cos(M_PI*rho/2.)*tgamma(-rho), 1./rho);
    double mu = -dt*p;

    int N = initialConditions.size();

    // Initialize

    Vector theta(N);
    for (int i = 0; i < N; ++i)
        theta(i) = initialConditions(i);

    Vector r(timeSteps+1);
    r(0) = order_param(theta);

    // Time stepping

    double U, V, W, xi, drift;
    double s = 0.0;
    long k = 0;

    do 
    {
        s += dt;
        k++;

        for (int i = 0; i < N; i++)
        {
            // simulate a tempered stable random variable

            do 
            {
                U = runif(rng);
                V = runif(rng);
                do {
                    W = rexp(rng);
                } while (W == 0.);
                xi = mu + sigma*cf*pow(A(M_PI*U,rho)/pow(W,1.-rho),1./rho);
            } while (V > exp(-alpha*xi));

            // calculate the drift

            drift = intrinsicFreqs(i);
            for (int j = 0; j < N; j++)
                drift -= globalCoupling/N*adjacencyMatrix(i,j)*sin(theta(i) - theta(j));

            // increment process

            theta(i) += drift*dt + noiseScale*xi;

            theta(i) = mod2pi(theta(i));
        }

        r(k) = order_param(theta);

#ifdef DEBUG
        printf("%6.3f", s);
        printf(" %6.3f", order_param(theta));
        for (int i = 0; i < N; i++) {
            printf(" %6.3f", theta(i));
        }
        printf("\n");
#endif

    } while (s < t);

    //std::cout << r << std::endl;
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
void paths(const double globalCoupling,
        const double noiseScale,
        Eigen::MatrixBase<DerivedA>& initialConditions,
        Eigen::MatrixBase<DerivedB>& intrinsicFreqs,
        Eigen::MatrixBase<DerivedC>& adjacencyMatrix,
        const double rho,
        const double c,
        const double alpha,
        const double p,
        const double t,
        const int timeSteps,
        const int seed,
        const int npaths) {

    std::mt19937 rng(seed);
    std::uniform_real_distribution<> runif(0,1);
    std::exponential_distribution<> rexp(1);

    double dt = t/timeSteps;

    double gamma = c/(1-rho) - p;
    double cf = pow(cos(M_PI_2*rho),-1./rho);
    double sigma = pow(-dt*c*cos(M_PI*rho/2.)*tgamma(-rho), 1./rho);
    double mu = -dt*p;

    int N = initialConditions.size();

    Vector avg_order(timeSteps+1);

    for (int j=0; j < npaths; ++j) {
        Vector theta(N);
        for (int i = 0; i < N; ++i)
            theta(i) = initialConditions(i);

        avg_order(0) += order_param(theta);

        double U, V, W, xi, drift;
        double s = 0.0;
        long k = 0;

        do {
            s += dt;
            k++;

            for (int i = 0; i < N; i++) {
                // simulate a tempered stable random variable
                do {
                    U = runif(rng);
                    V = runif(rng);
                    do {
                        W = rexp(rng);
                    } while (W == 0.);
                    xi = mu + sigma*cf*pow(A(M_PI*U,rho)/pow(W,1.-rho),1./rho);
                } while (V > exp(-alpha*xi));

                // calculate the drift
                drift = intrinsicFreqs(i);
                for (int j = 0; j < N; j++)
                    drift -= globalCoupling/N*adjacencyMatrix(i,j)*sin(theta(i) - theta(j));

                // increment process
                theta(i) += drift*dt + noiseScale*xi;
                theta(i) = mod2pi(theta(i));
            }

            avg_order(k) += order_param(theta);

        } while (s < t);
    }

    avg_order /= npaths;

    for (int i = 0; i < avg_order.size(); ++i) {
        printf("%6.3f %6.3f\n", i*dt, avg_order(i));
    }
}


int main(int argc, char const *argv[]) {
    double rho = 0.99; // stability index
    double c = 0.01;
    double alpha = 1.0;
    double xi = 0.2;
    double p = (1+xi)*(-c*rho*tgamma(-rho)*pow(alpha,rho-1));

    //printf("rho: %f c: %f alpha: %f xi: %f p: %f\n", rho, c, alpha, xi, p);

    double t = 30.0;
    int timeSteps = 5000;
    int seed = 1234;
    int npaths = 1000;

    double globalCoupling = 0.8;
    double noiseScale = 0.05;

    Vector initialConditions(3);
    initialConditions << 0.3, 6.0, 3.14;

    Vector intrinsicFreqs(3);
    intrinsicFreqs << 0.2, 0.2, 0.2;

    Matrix adjacencyMatrix(3,3);
    adjacencyMatrix <<
        0, 0, 1,
        0, 1, 0,
        0, 1, 1;

    path(globalCoupling, noiseScale, initialConditions, intrinsicFreqs,
            adjacencyMatrix, rho, c, alpha, p, t, timeSteps, seed);

    //paths(globalCoupling, noiseScale, initialConditions, intrinsicFreqs,
    //        adjacencyMatrix, rho, c, alpha, p, t, timeSteps, seed, npaths);

    return 0;
}
