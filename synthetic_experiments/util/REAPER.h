#ifndef REAPER_H
#define REAPER_H

#include <armadillo>

typedef struct{
    mat* X;
    vec* beta;
    int d;
} input_WLS;

typedef struct {
    mat P_star;
} output_WLS;

typedef struct{
    mat* X;
    int d;
} input_IRLS;

typedef struct {
    mat B_star;
} output_IRLS;

void WLS_solver(const input_WLS& in, output_WLS& out) {
    const int D = (*in.X).n_rows;
    const mat C = (*in.X) * diagmat(*in.beta) * (*in.X).t();
    vec v = zeros<vec>(D);

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, C);

    vec lambda = real(eigval);
    mat U = real(eigvec);
    uvec indices = sort_index(lambda, "descend");
    lambda = sort(lambda, "descend");
    U = U.cols(indices);

    double theta;
    if (lambda(in.d) == 0) {
        v.rows(0, in.d-1) = ones<vec>(in.d);
    } else {
        for (int i = in.d; i < D; ++i) {
            theta = (i-in.d+1) / sum(1/lambda.rows(0, i));
            if (i < D-1 && lambda(i) > theta && theta >= lambda(i+1)) {
                break;
            }
        }
        indices = find(lambda>theta);
        v(indices) = 1 - theta/lambda(indices);
    }

    out.P_star = U * diagmat(v) * U.t();
}


void IRLS_REAPER_solver(const input_IRLS& in, output_IRLS& out) {
    const int D = (*in.X).n_rows;
    const int L = (*in.X).n_cols;
    const double delta = 1e-10;
    const double epsilon = 1e-15;
    const int maxiter = 200;
    double alpha_old = datum::inf;
    vec beta = ones<vec>(L);

    input_WLS inWLS;
    inWLS.X = in.X;
    inWLS.d = in.d;
    output_WLS outWLS;
    mat temp1;
    vec temp2;
    double alpha_new;
    int i = 0;
    while(true) {
        inWLS.beta = &beta;
        WLS_solver(inWLS, outWLS);
        temp1 = *in.X - outWLS.P_star * (*in.X);
        temp1 = sqrt(sum(temp1 % temp1, 0));
        temp2 = temp1.t();
        beta = 1 / max(delta*ones<vec>(L), temp2);
        alpha_new = sum(temp2);
        if (alpha_new >= alpha_old - epsilon || i++ > maxiter) {
            break;
        }
        alpha_old = alpha_new;
    }

    mat U, V;
    vec s;
    svd(U,s,V, outWLS.P_star);

    // out.B_star = U.cols(in.d, D-1);
    out.B_star = U.col(D-1);
}


#endif
