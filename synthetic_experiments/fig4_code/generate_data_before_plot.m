clc
clear

addpath ../util/

D = 30;
d = 29;

% set sigma = 0 to plot a; sigma = 0.1 to plot b
sigma = 0;
for t = 1:5
    N_val{t} = [];
    M_val{t} = [];
    cos_phi{t} = [];
    for N = 10:10:500
        barX = [randn(d,N); zeros(D-d, N)]/sqrt(d);
        for M = 10:50:3000
            O = randn(D, M)/sqrt(D); O = normc(O);
            barE = sigma * randn(D, N) / sqrt(D);
            v_norm = vecnorm(barX+barE);
            m_norm = repmat(v_norm, D, 1);
            X = barX ./ m_norm;
            noise = barE ./ m_norm;

            [b,~] = DPCP_PSGM(X,O,noise);
            cos_phi{t} = [cos_phi{t}, norm(b(1:d))];
            N_val{t} = [N_val{t}, N];
            M_val{t} = [M_val{t}, M];
        end
    end
end


N = N_val{1};
M = M_val{1};
cos_phi = (cos_phi{1}+cos_phi{2}+cos_phi{3}+cos_phi{4}+cos_phi{5})/5;

save('M.txt', 'M', '-ascii')
save('N.txt', 'N', '-ascii')
save('cos_phi.txt', 'cos_phi', '-ascii')

