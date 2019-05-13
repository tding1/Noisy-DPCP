#ifndef PSGD_H
#define PSGD_H

#include <armadillo>

using namespace arma;

typedef struct {
    mat* X;
    bool minimize;
} input_PSGD;

typedef struct {
    mat* X;
    vec* b;
    bool minimize;
} input_obj;

typedef struct {
    vec b_star;
    double optimum;
} output_PSGD;

double obj(const input_obj& in) {
    double temp = norm((*in.X).t() * (*in.b), 1);
    return (in.minimize == true) ? temp : -temp;
}

void DPCP_PSGM(const input_PSGD& in, output_PSGD& out) {

    input_obj in_obj;
    const double mu_min      = 1e-15;
    const double mu_0        = 1e-2;
    const int maxiter        = 200;
    const double alpha       = 1e-3;
    const double beta        = 0.5;

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, (*in.X) * (*in.X).t());
    uword ind = (in.minimize == true) ?  eigval.index_min() : eigval.index_max();
    vec b = real(eigvec.col(ind));
    b /= norm(b);
    
    in_obj.X        = in.X;
    in_obj.b        = &b;
    in_obj.minimize = in.minimize;
    double obj_old  = obj(in_obj);
    double mu       = mu_0;
    double grad_norm_square;
    vec b_next, grad;
    int i = 0;
    while (mu > mu_min && i <= maxiter) {
        i++;
        grad = (*in.X) * sign((*in.X).t() * b);
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
