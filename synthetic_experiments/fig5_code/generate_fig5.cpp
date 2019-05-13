#include <iostream>
#include <cstdlib>
#include "../util/PSGD.h"
#include "../util/piecewise_PSGD.h"
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

int main(int argc, const char** argv) {
    wall_clock timer;
    timer.tic();
    
    // arma_rng::set_seed_random();
    arma_rng::set_seed(1234);

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
    const double sigma_limit = cmdInput[5];//0.1;

    const int M = ceil(r * N / (1 - r));
    mat init_X = 1/sqrt(d)*join_vert(randn<mat>(d, N), zeros<mat>(D-d, N));
    mat O = 1/sqrt(D)*randn<mat>(D, M); O = normalise(O);

    // simulate eta_O
    double eta_O = -datum::inf;
    double cp, sp, temp;
    vec sample_b, proj_b, xi_b, s, n, diff;
    for (int i = 0; i < 1e5; ++i) {
        sample_b = randn<vec>(D);
        sample_b /= norm(sample_b);
        proj_b = join_vert(sample_b.rows(0, d-1), zeros<vec>(D-d));

        cp = norm(proj_b);
        sp = sqrt(1 - pow(cp, 2));
        s = proj_b / cp;
        diff = sample_b - proj_b;
        n = diff / norm(diff);
        xi_b = -sp * s + cp * n;

        temp = abs(sum((sign(O.t() * sample_b)) % (O.t() * xi_b)));
        if (temp > eta_O) {
            eta_O = temp;
        }
    }

    const double sigma   = sigma_limit;;
    double chatxmin      = 0;
    double chatxmax      = 0;
    double chatemax      = 0;

    input_PSGD in;
    output_PSGD out;

    mat Y = O;
    in.X = &Y;
    in.minimize = false;
    DPCP_PSGM(in, out);
    double comax = out.optimum;
    in.minimize = true;
    DPCP_PSGM(in, out);
    double comin = out.optimum;

    mat noise = sigma * randn<mat>(D, N);
    mat sum_x_noise = init_X + noise;
    mat norm_col = ones<mat>(D, N);
    for (int i = 0; i < D; i++)
        for (int j = 0; j < N; j++)
            norm_col(i,j) = norm(sum_x_noise.cols(j,j));

    mat X = init_X / norm_col;
    noise = noise / norm_col;
    mat inlier_noise = join_vert(noise.rows(0, d-1), zeros<mat>(D-d, N));
    mat outlier_noise = noise - inlier_noise;

    mat Xhat = X + inlier_noise;
    mat Ehat = outlier_noise;

    Y = Xhat.rows(0, d-1);
    in.X = &Y;
    in.minimize = true;
    DPCP_PSGM(in, out);
    chatxmin = out.optimum;

    Y = Xhat.rows(0, d-1);
    in.X = &Y;
    in.minimize = false;
    DPCP_PSGM(in, out);
    chatxmax = out.optimum;

    Y = Ehat.rows(d, D-1);
    in.X = &Y;
    in.minimize = false;
    DPCP_PSGM(in, out);
    chatemax = out.optimum;
    
    double xi_hate = sqrt(N) * norm(Ehat);

    double m = -1;
    for (int k = 1; k <= 2e4; k++) {
        vec s = randn<vec>(D);
        s.rows(d,D-1) = zeros<vec>(D-d);
        s /= norm(s);
        
        vec n = randn<vec>(D);
        n.rows(0,d-1) = zeros<vec>(d);
        n /= norm(n);

        double tmp = norm((eye(D,D)-s*s.t()) * Xhat * sign(Xhat.t() * s + randi<int>(distr_param(0, 10000)) * Ehat.t() * n));
        if (tmp > m) {
            m = tmp;
        }
    }

    vec temp2 = {chatxmax, comax, xi_hate};
    double mu_prime = 1.0 / (16 * max(temp2));
    double theta_prime = atan(1.2*(chatemax+xi_hate)/(chatxmin-m-eta_O-xi_hate));
    double theta_bound = atan((chatxmin-chatemax-xi_hate)/(m+eta_O+xi_hate));
    cout << "noise: " << sigma << endl;
    cout << "theta0 bound: " << theta_bound * 180 / datum::pi << endl;
    cout << "mu' = " << mu_prime << endl;
    cout << "theta' = " << theta_prime * 180 / datum::pi << endl;
    cout << "Inlier requirement holds ? " << (chatxmin > m + eta_O + chatemax * 2.5 + xi_hate * 3.5) << endl;

    mat Xtilde_noise = join_horiz(X + noise, O);

    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, Xtilde_noise * Xtilde_noise.t());
    uword ind = eigval.index_min();
    vec init_b = real(eigvec.col(ind));
    init_b /= norm(init_b);
    double theta0 = asin(norm(init_b.rows(0, d-1)));
    
    cout << "theta0 = " << theta0 * 180 / datum::pi << endl;

    int K0 = 5;
    double mu0 = mu_prime;
    if (theta0 >= theta_prime) {
        double t1 = chatxmin-m-eta_O-xi_hate-1/tan(theta0)*(chatemax+xi_hate);
        double t2 = chatxmin- tan(theta0) * (m+eta_O+xi_hate) - (chatemax+xi_hate);
        K0 = tan(theta0) / mu0 / ((t1<=t2) ? t1 : t2) + 1;
    }
    cout << "K0 = " << K0 << endl;

    double beta = 0.5;
    int K = 42 / (25*beta*mu_prime*(chatxmin-m-eta_O-xi_hate))*0.5+1;
    cout << "K = " << K << endl;

    input_PPSGD inP;
    output_PPSGD outP;
    inP.X = &Xtilde_noise;
    inP.minimize = true;
    inP.init_b = &init_b;
    inP.K0 = K0;
    inP.K = K;
    inP.mu0 = mu0;
    inP.beta = beta;
    inP.theta0 = theta0;
    inP.theta_prime = theta_prime;
    inP.d = d;
    inP.comin = comin;
    PPSGM(inP, outP);
    vec stepsize = outP.stepsize;
    vec tan_theta = outP.tan_theta;
    vec bound = outP.bound;

    stepsize.save("./files/stepsize.ty", raw_ascii);
    tan_theta.save("./files/tan_theta.ty", raw_ascii);
    bound.save("./files/bound.ty", raw_ascii);

    in.X = &Xtilde_noise;
    in.minimize = true;
    DPCP_PSGM(in, out);
    // cout << outP.optimum << " " << out.optimum << endl;

    cout << tan(asin(norm(out.b_star.rows(0, d-1)))) << endl;

    cout << endl;

    cout << "Data generated... Time elapsed " << timer.toc() << "s." << endl;


    
    cout << "Done!" << endl;

    return 0;
}
