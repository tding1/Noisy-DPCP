#ifndef GENERAL_PSGD_H
#define GENERAL_PSGD_H

#include <armadillo>

using namespace arma;

typedef struct {
    mat* A;
    mat* B;
    bool minimize;
} input_generalPSGD;

typedef struct {
    mat* A;
    mat* B;
    vec* b;
    bool minimize;
} input_obj_generalPSGD;

typedef struct {
    vec b_star;
    double optimum;
} output_generalPSGD;

double obj(const input_obj_generalPSGD& in) {
    double temp = norm((*in.A).t() * (*in.b), 1) - norm((*in.B).t() * (*in.b), 1);
    return (in.minimize == true) ? temp : -temp;
}

void DPCP_generalPSGDone(const input_generalPSGD& in, output_generalPSGD& out) {

    input_obj_generalPSGD in_obj;
    const double mu_min      = 1e-15;
    const double mu_0        = 1e-2;
    const double maxiter     = 200;
    const double alpha       = 1e-3;
    const double beta        = 0.5;

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, (*in.A) * (*in.A).t());
    uword ind = (in.minimize == true) ?  eigval.index_min() : eigval.index_max();
    vec b = real(eigvec.col(ind));
    b /= norm(b);
    
    in_obj.A        = in.A;
    in_obj.B        = in.B;
    in_obj.b        = &b;
    in_obj.minimize = in.minimize;
    double obj_old  = obj(in_obj);
    double mu       = mu_0;
    double grad_norm_square;
    vec b_next, grad;
    int i = 0;
    while (mu > mu_min && i <= maxiter) {
        i++;
        grad = (*in.A) * sign((*in.A).t() * b) - (*in.B) * sign((*in.B).t() * b);
        if (!in.minimize) {
            grad = -grad;
        }
        grad_norm_square = pow(norm(grad), 2);

        b_next = b - mu * grad;
        b_next /= norm(b_next);
        while (true) {
            if (mu <= mu_min)
                break;
            in_obj.b = &b_next;
            if (obj(in_obj) <= obj_old - alpha * mu * grad_norm_square)
                break;
            mu *= beta;
            b_next = b - mu * grad;
            b_next /= norm(b_next);
        }
        b = b_next;
        in_obj.b = &b;
        obj_old = obj(in_obj);
    }

    out.b_star = b;
    out.optimum = (in.minimize == true) ? obj_old : -obj_old;
}

#endif
