#ifndef GGD_H
#define GGD_H

#include <armadillo>

typedef struct{
    mat* X;
    int d;
} input_GGD;

typedef struct {
    mat V;
} output_GGD;

typedef struct {
    mat* X;
    mat* V;
    int D;
} input_obj_GGD;

double obj(const input_obj_GGD& in) {
    vec tmp = sum(square((eye<mat>(in.D, in.D) - (*in.V) * (*in.V).t()) * (*in.X)), 0).t();
    return sum(sqrt(tmp));
}

void GGD_solver(const input_GGD& in, output_GGD& out) {
    wall_clock timer;
    input_obj_GGD in_obj;
    const double mu_min      = 1e-15;
    const double mu_0        = 1e-2;
    const int maxiter        = 200;
    const double alpha       = 1e-3;
    const double beta        = 0.5;
    const int D              = (*in.X).n_rows;
    const int d              = in.d;

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, (*in.X) * (*in.X).t());

    vec lambda = real(eigval);
    mat V = real(eigvec);
    uvec indices = sort_index(lambda, "descend");
    V = V.cols(indices);
    V = V.cols(0, d-1);

    in_obj.X        = in.X;
    in_obj.V        = &V;
    in_obj.D        = D;
    double obj_old  = obj(in_obj);
    double mu       = mu_0;
    double grad_norm_square;
    mat V_next, QV, partial, grad, Y, T, U, W;
    vec tmp, tmp2, s;
    uvec ind;
    int i = 0;
    // double ttt = 0;
    while (mu > mu_min && i <= maxiter) {
        i++;
        QV = eye<mat>(D, D) - V * V.t();
        tmp = sqrt(sum(square(QV * (*in.X)), 0).t());
        ind = find(tmp > 0);
        Y = (*in.X).cols(ind);
        tmp2 = tmp.rows(ind);
    
        // timer.tic();
        T = sqrt(repmat(tmp2, 1, D));
        partial = -(Y / T.t()) * (Y.t() / T) * V;
        // ttt += timer.toc();

        grad = QV * partial;
        grad_norm_square = pow(norm(grad, "fro"), 2);
        svd_econ(U, s, W, -grad);

        do {
            if (mu <= mu_min)
                break;
            V_next = V*W*diagmat(cos(s*mu))*W.t() + U*diagmat(sin(s*mu))*W.t();
            in_obj.V = &V_next;
            if (obj(in_obj) <= obj_old - alpha * mu * grad_norm_square)
                break;
            mu *= beta;
        } while (true);

        V = V_next;
        in_obj.V = &V;
        obj_old = obj(in_obj);
    }

    out.V = V;
}


#endif
