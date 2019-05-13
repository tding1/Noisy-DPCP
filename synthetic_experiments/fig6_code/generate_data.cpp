#include <iostream>
#include <cstdlib>
#include "../util/PSGD.h"
#include <armadillo>

using namespace arma;
using std::cout;
using std::endl;

void parse_Cmdline(const char** argv, double* input) {
    for (int i = 1; i < 7; i += 2) {
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
        }
    }
}

int main(int argc, const char** argv) {
    wall_clock timer;
    timer.tic();
    
    arma_rng::set_seed_random();
    // arma_rng::set_seed(1234);

    if (argc != 7) {
        cout << "Not enough input!" << endl;
        return 1;
    }

    double cmdInput[3];
    parse_Cmdline(argv, cmdInput);

    const int D              = (int)(cmdInput[0]); //30;
    const int d              = (int)(cmdInput[1]);//29;
    const int N              = (int)(cmdInput[2]);//1500;

    double ratio[8] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7};
    mat success = zeros<mat>(8, 8);
    mat success_repear = zeros<mat>(8, 8);
    for (int rr = 0; rr < 8; rr++) {
        double r = ratio[rr];
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

        
        const vec sigma      = {0, 0.01,0.02,0.03,0.04,0.05,0.06, 0.07};
        const int num_seg    = 8;
        vec chatxmin         = zeros<vec>(num_seg);
        vec chatemax         = zeros<vec>(num_seg);
        vec comin            = zeros<vec>(num_seg);
        vec comax            = zeros<vec>(num_seg);
        vec alpha            = zeros<vec>(num_seg);
        vec beta             = zeros<vec>(num_seg);
        vec cos_phi          = zeros<vec>(num_seg);
        vec t1               = zeros<vec>(num_seg);
        vec t2               = zeros<vec>(num_seg);
        vec approx_roots     = zeros<vec>(num_seg);
        vec etaO             = eta_O * ones<vec>(num_seg);
        mat noise, inlier_noise, outlier_noise, Xhat, Ehat, Xtilde_noise, Y, Z;
        vec coeff, rts, hat_b;
        double alp2, be2, F;
        input_PSGD in;
        output_PSGD out;
        int num_test = 10;
        for (int tt = 0; tt < num_test; tt++) {
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

                Y = O;
                in.X = &Y;
                in.minimize = true;
                DPCP_PSGM(in, out);
                comin[i] = out.optimum;

                Y = O;
                in.X = &Y;
                in.minimize = false;
                DPCP_PSGM(in, out);
                comax[i] = out.optimum;

                beta[i] = chatemax[i] / chatxmin[i];
                alpha[i] = (eta_O + D) / chatxmin[i];

                alp2 = pow(alpha[i], 2);
                be2 = pow(beta[i], 2);
                coeff = {1, 0, alp2 - 1, 4*alpha[i]*beta[i], 4*be2};
                rts = real(roots(coeff));
                rts = sort(rts);
                t1[i] = rts[2];
                t2[i] = rts[3];

                F = chatxmin[i] / (4*sqrt(d)) - norm(O) * norm(normalise(join_vert(zeros<mat>(d, M), O.rows(d, D-1))));
                temp = 0;
                for (int j = 0; j < N; ++j) {
                    temp += norm(Ehat.col(j));
                }

                if (F-temp > 0) {
                    success_repear(rr, i)++;
                }

                if (((comax[i] - comin[i])/chatxmin[i]+2*beta[i] < t2[i])
                    && (pow(sqrt(alp2+8) - 3*alpha[i], 1.5) * pow(sqrt(alp2+8) + alpha[i], 0.5) >=  32*beta[i] )) {
                    success(rr, i)++;
                }
            }
        }
    }

    cout << "Data generated... Time elapsed " << timer.toc() << "s." << endl;
    cout << "Done!" << endl;
    success = success/10;
    cout << success << endl;
    success.save("success.ty", raw_ascii);

    success_repear = success_repear/10;
    cout << success_repear << endl;
    success_repear.save("success_reaper.ty", raw_ascii);

    return 0;
}
