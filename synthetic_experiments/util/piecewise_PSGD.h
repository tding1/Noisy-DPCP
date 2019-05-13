#ifndef PIECEWISEPSGD_H
#define PIECEWISEPSGD_H

#include <armadillo>
#include "PSGD.h"

using std::cout;
using std::endl;

using namespace arma;

typedef struct {
    mat* X;
    vec* init_b;
    int K0;
    int K;
    double mu0;
    double beta;
    double theta_prime;
    double theta0;
    double comin;
    int d;
    bool minimize;
} input_PPSGD;

typedef struct {
    vec b_star;
    vec bound;
    vec tan_theta;
    vec stepsize;
    double optimum;
} output_PPSGD;

void PPSGM(const input_PPSGD& in, output_PPSGD& out) {

    const double mu_min      = 1e-30;
    const double mu_0        = in.mu0;
    const double beta        = in.beta;
    const int K0             = in.K0;
    const int K              = in.K;
    const int maxiter        = 5000;
    const double theta0      = in.theta0;
    const double theta_prime = in.theta_prime;

    vec bound                = zeros<vec>(maxiter);
    vec tan_theta            = zeros<vec>(maxiter);
    vec stepsize             = zeros<vec>(maxiter);

    vec b = *in.init_b;

    double t1 = tan(theta0);
    double t2 = 7.0/25 + 1.5 * sin(theta_prime) * (1-mu_0*in.comin);
    double init_bound = (t1 <= t2) ? t2 : t1;
    
    input_obj in_obj;
    in_obj.X        = in.X;
    in_obj.b        = in.init_b;
    in_obj.minimize = in.minimize;
    double obj_old  = obj(in_obj);
    vec b_next, grad;
    double mu = 1;
    int i = 0;
    while (mu > mu_min && i <= maxiter) {
        i++;
        grad = (*in.X) * sign((*in.X).t() * b);
        if (!in.minimize) {
            grad = -grad;
        }


        if (i < K0){
            mu = mu_0;
            bound[i-1] = init_bound;
        } else {
            mu = mu_0 * pow(beta, (i-K0)/K + 1);
            bound[i-1] = 7.0/25 * pow(beta, (i-K0)/K) + 1.5 * sin(theta_prime) * (1-mu*in.comin);
        }

        stepsize[i-1] = mu;

        b_next = b - mu * grad;
        b_next /= norm(b_next);
        tan_theta[i-1] = tan(asin(norm(b_next.rows(0, in.d-1))));
    
        b = b_next;
        in_obj.b = &b;
        obj_old = obj(in_obj);
    }

    out.b_star = b;
    out.optimum = (in.minimize == true) ? obj_old : -obj_old;
    out.stepsize = stepsize;
    out.tan_theta = tan_theta;
    out.bound = bound;
}

#endif
