#include "variates.h"

int main(int argc, char const *argv[]) {
    int seed      = argc > 1 ? atol(argv[1]) : 1234;
    int npaths    = argc > 2 ? atoi(argv[2]) : 1000;
    int nsteps    = argc > 3 ? atoi(argv[3]) : 5000;
    double alpha  = argc > 4 ? atof(argv[4]) : 1.7; // stability
    double lambda = argc > 5 ? atof(argv[5]) : 1.0; // tempering
    double sigma  = argc > 6 ? atof(argv[6]) : 1.0; // diffusivity
    double K      = argc > 7 ? atof(argv[7]) : 0.8; // global coupling
    double max_t  = argc > 8 ? atof(argv[8]) : 30.0; // maximum time
    
    double dt = max_t/nsteps;
    double a = 0.5*alpha*pow(sigma,2.0)/(tgamma(1-alpha)*cos(M_PI*alpha/2));
    double b = lambda;

    TemperedStableDistribution rtstable(alpha, dt*a, b, 1.1);
    NormalDistribution rnorm(0, dt*pow(sigma,2));
    mt19937 rng(seed);
    
    for (int i = 0; i < 5000; i++) {
        if (alpha > 1.999) {
            printf("%.4f\n", rnorm(rng));
        } else {
            printf("%.4f\n", rtstable(rng));
        }
    }
   
    return 0;
}
