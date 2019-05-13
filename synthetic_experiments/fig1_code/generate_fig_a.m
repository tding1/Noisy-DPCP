clc
clear

addpath ../util/

rng(2018);

D = 30;
N = 500;
r = 0.7;
M = ceil(r * N / (1 - r));
d = 29;
sigma = .05;
segnum = 100;

barX = [randn(d,N); zeros(D-d, N)]/sqrt(d);
O = randn(D, M)/sqrt(D); O = normc(O);

barE = sigma * randn(D, N) / sqrt(D);
v_norm = vecnorm(barX+barE);
m_norm = repmat(v_norm, D, 1);
X = barX ./ m_norm;
noise = barE ./ m_norm;

[b,~] = DPCP_PSGM(X,O,noise);
cos_phi_DPCP = norm(b(1:d));

tau = linspace(0.001,0.15,segnum);
cos_denoised = zeros(segnum,1);
for i = 1:size(tau,2)
    [b, ~, ~] = denoisedDPCP(X,O,noise, tau(i));
    cos_denoised(i) = norm(b(1:d));
end

plot(tau, cos_denoised, tau, cos_phi_DPCP * ones(100,1), 'r:', 'linewidth', 4)
legend({'denoised version (2)', 'formulation (3)'}, 'FontSize', 30)
xlabel('\tau','FontSize', 32)
ylabel('$\sin(\theta_*)$','FontSize', 32,'FontName','Times New Roman','Interpreter','LaTex')


set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , 35             , ...
    'FontName'  , 'Times New Roman');
