#include <iostream>
#include <cstdlib>
#include "../util/PSGD.h"
#include "../util/REAPER.h"
#include "../util/GGD.h"
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
    const mat X = join_vert(randn<mat>(d, N), zeros<mat>(D-d, N));
    const mat O = randn<mat>(D, M);

    const vec sigma      = linspace<vec>(0, sigma_limit, num_seg);
    vec cos_phi_PSGD     = zeros<vec>(num_seg);
    vec cos_phi_REAPER   = zeros<vec>(num_seg);
    vec cos_phi_GGD      = zeros<vec>(num_seg);
    mat noise, Xtilde_noise, Y;

    input_PSGD in_PSGD;
    output_PSGD out_PSGD;
    input_IRLS in_IRLS;
    output_IRLS out_IRLS;
    input_GGD in_GGD;
    output_GGD out_GGD;
    for (int i = 0; i < num_seg; ++i) {
        noise = sigma[i] * randn<mat>(D, N);
        Xtilde_noise = join_horiz(X + noise, O);

        // timer.tic();
        in_PSGD.X = &Xtilde_noise;
        in_PSGD.minimize = true;
        DPCP_PSGM(in_PSGD, out_PSGD);
        cos_phi_PSGD[i] = norm(out_PSGD.b_star.rows(0, d-1));
        // cout << timer.toc() << endl;

        // timer.tic();
        in_IRLS.X = &Xtilde_noise;
        in_IRLS.d = d;
        IRLS_REAPER_solver(in_IRLS, out_IRLS);
        cos_phi_REAPER[i] = norm(out_IRLS.B_star.rows(0, d-1));
        // cout << timer.toc() << endl;

        // timer.tic();
        in_GGD.X = &Xtilde_noise;
        in_GGD.d = d;
        GGD_solver(in_GGD, out_GGD);
        Y = null((out_GGD.V).t());
        cos_phi_GGD[i] = norm(Y.rows(0, d-1));
        // cout << timer.toc() << endl;

    }

    cout << "Data generated... Time elapsed " << timer.toc() << "s." << endl;

    cos_phi_REAPER.save("./files/cos_phi_REAPER.ty", raw_ascii);
    cos_phi_PSGD.save("./files/cos_phi_PSGD.ty", raw_ascii);
    cos_phi_GGD.save("./files/cos_phi_GGD.ty", raw_ascii);

    cout << "Done!" << endl;

    return 0;
}
