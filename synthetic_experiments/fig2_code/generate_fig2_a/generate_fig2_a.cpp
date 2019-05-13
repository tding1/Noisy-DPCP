#include <iostream>
#include <cstdlib>
#include "../../util/PSGD.h"
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
    const int num_seg        = (int)(cmdInput[4]);//200;
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

    const vec sigma      = linspace<vec>(0, sigma_limit, num_seg);
    vec chatxmin         = zeros<vec>(num_seg);
    vec chatemax         = zeros<vec>(num_seg);
    vec alpha            = zeros<vec>(num_seg);
    vec beta             = zeros<vec>(num_seg);
    vec etaO             = eta_O * ones<vec>(num_seg);

    input_PSGD in;
    output_PSGD out;
    for (int i = 0; i < num_seg; ++i) {
        mat noise = 1/sqrt(D)*sigma[i] * randn<mat>(D, N);

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

        mat Y = Xhat.rows(0, d-1);
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
    }

    cout << "Data generated... Time elapsed " << timer.toc() << "s." << endl;

    chatxmin.save("./files/chatxmin.ty", raw_ascii);
    chatemax.save("./files/chatemax.ty", raw_ascii);
    etaO.save("./files/etaO.ty", raw_ascii);
    alpha.save("./files/alpha.ty", raw_ascii);
    beta.save("./files/beta.ty", raw_ascii);
    
    cout << "Done!" << endl;

    return 0;
}