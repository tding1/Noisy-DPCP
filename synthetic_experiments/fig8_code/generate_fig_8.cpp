#include <iostream>
#include <cstdlib>
#include "../util/PSGD.h"
#include <armadillo>

using namespace arma;
using std::cout;
using std::endl;

void parse_Cmdline(const char** argv, double* input) {
    for (int i = 1; i < 12; i += 2) {
        switch(argv[i][1]) {
            case 'D':
                input[0] = atof(argv[i+1]);
                break;
            case 'd':
                input[1] = atof(argv[i+1]);
                break;
            case 'N':
                input[2] = atof(argv[i+1]);
                break;
            case 'r':
                input[3] = atof(argv[i+1]);
                break;
            case 'n':
                input[4] = atof(argv[i+1]);
                break;
            case 's':
                input[5] = atof(argv[i+1]);
                break;
        }
    }
}

void findGammaAndEta(double* gamma, double* eta, const int& d, const mat& Xtilde_noise) {
    const double increment = 0.05;
    const int D = Xtilde_noise.n_rows;
    const int L = Xtilde_noise.n_cols;
    mat X1, X0, T;
    vec v;
    double P, A, S;
    cx_vec eigval;
    cx_mat eigvec;

    do {
        X1.clear();
        X0.clear();
        T = zeros<mat>(d, d);

        *gamma -= increment; 
        *eta += increment;
        for (int j = 0; j < L; j++) {
            v = Xtilde_noise.col(j);
            if (norm(v.rows(d, D-1)) / norm(v) <= sin(*eta))
                X1 = join_horiz(X1, v);
            else
                X0 = join_horiz(X0, v);
        }

        for (uword j = 0; j < X1.n_cols; j++) {
            v = X1.col(j);
            v = v.rows(0, d-1);
            T += v * v.t() / norm(v);
        }
        eig_gen(eigval, eigvec, T);
        P = min(real(eigval));

        A = sqrt(X0.n_cols) * norm(X0);

        S = cos(*gamma + *eta) * cos(*eta) * P - A;
        cout << "    S = " << S << " P = " << P << " A = " << A << endl;
        if (S > 0) break;

    } while ((*gamma > 0) && (*eta < datum::pi / 2) && (*eta < *gamma));
   
    if (S <= 0 || *eta > *gamma) {
        *gamma = -1;
        *eta = -1;
    }

}

int main(int argc, const char** argv) {
    wall_clock timer;
    timer.tic();
    
    arma_rng::set_seed_random();
    // arma_rng::set_seed(1234);

    if (argc != 13) {
        cout << "Not enough input!" << endl;
        return 1;
    }

    double cmdInput[6];
    parse_Cmdline(argv, cmdInput);

    const int D              = (int)(cmdInput[0]); //30;
    const int d              = (int)(cmdInput[1]);//29;
    const int N              = (int)(cmdInput[2]);//1500;
    const double r           = cmdInput[3];//0.7;
    const int num_seg        = (int)(cmdInput[4]);//200;
    const double sigma_limit = cmdInput[5];//0.1;

    const int M = ceil(r * N / (1 - r));
    const int L = M + N;
    mat init_X = 1/sqrt(d)*join_vert(randn<mat>(d, N), zeros<mat>(D-d, N));
    mat O = 1/sqrt(D)*randn<mat>(D, M); O = normalise(O);

    // simulate eta_O
    double eta_O = -datum::inf;
    double temp;
    vec sample_b, proj_b, xi_b, s, n, diff;
    for (int i = 0; i < 1e5; ++i) {
        sample_b = randn<vec>(D);
        sample_b /= norm(sample_b);

        mat G = null(sample_b.t());
        xi_b = G.col(randi<int>(distr_param(0,G.n_cols-1)));
        mat prod = xi_b.t() * O * sign(O.t() * sample_b);
        temp = abs(prod(0,0));

        if (temp > eta_O) {
            eta_O = temp;
        }
    }
    
    const vec sigma      = linspace<vec>(0, sigma_limit, num_seg);
    vec chatxmin         = zeros<vec>(num_seg);
    vec chatemax         = zeros<vec>(num_seg);
    vec alpha            = zeros<vec>(num_seg);
    vec beta             = zeros<vec>(num_seg);
    vec t1               = zeros<vec>(num_seg);
    vec t2               = zeros<vec>(num_seg);
    vec angle_eta        = zeros<vec>(num_seg);
    vec angle_gamma      = zeros<vec>(num_seg);
    vec theta2           = zeros<vec>(num_seg);
    vec hat_theta2       = zeros<vec>(num_seg);

    mat noise, inlier_noise, outlier_noise, Xtilde_noise, Xhat, Ehat, Y;
    mat X1, X0, T, G;
    vec rts, v;
    input_PSGD in;
    output_PSGD out;
    double alp2, be2, GGD_eta_O, P;
    cx_vec eigval; cx_mat eigvec;
    for (int i = 0; i < num_seg; ++i) {
        noise = 1/sqrt(D)*sigma[i] * randn<mat>(D, N);

        mat sum_x_noise = init_X + noise;
        mat norm_col = ones<mat>(D, N);
        for (int i = 0; i < D; i++)
            for (int j = 0; j < N; j++)
                norm_col(i,j) = norm(sum_x_noise.cols(j,j));

        mat X = init_X / norm_col;
        noise = noise / norm_col;
        
        inlier_noise = join_vert(noise.rows(0, d-1), zeros<mat>(D-d, N));
        outlier_noise = noise - inlier_noise;

        Xhat = X + inlier_noise;
        Ehat = outlier_noise;

        Y = Xhat.rows(0, d-1);
        in.X = &Y;
        in.minimize = true;
        DPCP_PSGM(in, out);
        chatxmin[i] = out.optimum;

        Y = Ehat.rows(d, D-1);
        in.X = &Y;
        in.minimize = false;
        DPCP_PSGM(in, out);
        chatemax[i] = out.optimum;

        beta[i] = chatemax[i] / chatxmin[i];
        alpha[i] = eta_O / chatxmin[i];

        alp2 = pow(alpha[i], 2);
        be2 = pow(beta[i], 2);
        vec coeff = {1, 0, alp2 - 1, 4*alpha[i]*beta[i], 4*be2};
        rts = real(roots(coeff));
        rts = sort(rts);
        t1[i] = rts[2];
        t2[i] = rts[3];

        angle_eta[i] = asin(t1[i]);
        theta2[i] = asin(t2[i]) * 180 / datum::pi;

        Xtilde_noise = join_horiz(X + noise, O);
        X1.clear(); X0.clear();
        for (int j = 0; j < L; j++) {
            v = Xtilde_noise.col(j);
            if (norm(v.rows(d, D-1)) / norm(v) <= sin(angle_eta[i]))
                X1 = join_horiz(X1, v);
            else
                X0 = join_horiz(X0, v);
        }

        GGD_eta_O = -datum::inf;
        for (int i = 0; i < 1e5; ++i) {
            sample_b = randn<vec>(D);
            sample_b /= norm(sample_b);

            G = null(sample_b.t());
            xi_b = G.col(randi<int>(distr_param(0,G.n_cols-1)));
            mat prod = xi_b.t() * X0 * sign(X0.t() * sample_b);
            temp = abs(prod(0,0));

            // temp = abs(sum((sign(X0.t() * sample_b)) % (X0.t() * xi_b)));
            if (temp > GGD_eta_O) {
                GGD_eta_O= temp;
            }
        }

        T = zeros<mat>(d, d);
        for (uword j = 0; j < X1.n_cols; j++) {
            v = X1.col(j);
            v = v.rows(0, d-1);
            temp = norm(v);
            if (temp < 1e-6)
                continue;
            T += v * v.t() / temp;
        }
        eig_gen(eigval, eigvec, T);
        P = min(real(eigval));
        
        temp = GGD_eta_O / (cos(angle_eta[i]) * P);
        hat_theta2[i] = (acos(temp) - angle_eta[i]) * 180 / datum::pi;

        cout << "sigma = " << sigma[i] << ", arcsin(t2) = " << theta2[i] << ", gamma = " << hat_theta2[i] << endl;
    }

    theta2.save("./files/theta2.ty", raw_ascii);
    hat_theta2.save("./files/hat_theta2.ty", raw_ascii);

    return 0;
}
