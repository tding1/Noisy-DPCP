clc
clear

addpath ../util/

rng(2018);

D = 30;
N = 500;
r = 0.7;
M = ceil(r * N / (1 - r));
d = 29;

barX = [randn(d,N); zeros(D-d, N)]/sqrt(d);
O = randn(D, M)/sqrt(D); O = normc(O);

num_seg = 100;
sigma = linspace(0, 0.2, num_seg);
tau = [0.005, 0.015, 0.03];
cos_phi_DPCP = zeros(num_seg, 1);
cos_phi_denoised = zeros(num_seg, size(tau,2));

for i = 1:num_seg
    barE = sigma(i) * randn(D, N) / sqrt(D);
    
    v_norm = vecnorm(barX+barE);
    m_norm = repmat(v_norm, D, 1);
    X = barX ./ m_norm;
    noise = barE ./ m_norm;
    
    [b,~] = DPCP_PSGM(X,O,noise);
    cos_phi_DPCP(i) = norm(b(1:d));
    
    for j = 1:size(tau,2)
        [b, ~, ~] = denoisedDPCP(X,O,noise, tau(j));
        cos_phi_denoised(i, j) = norm(b(1:d));
    end

end

plot(sigma, cos_phi_DPCP, 'linewidth',4, 'DisplayName','formulation (3)')
hold on
plot(sigma, cos_phi_denoised(:, 1), ':', 'linewidth',4, 'DisplayName', strcat('\tau = ', num2str(tau(1))))
plot(sigma, cos_phi_denoised(:, 2), '--', 'linewidth',4, 'DisplayName', strcat('\tau = ', num2str(tau(2))))
plot(sigma, cos_phi_denoised(:, 3), '-.', 'linewidth',4, 'DisplayName', strcat('\tau = ', num2str(tau(3))))

lgd = legend;
lgd.FontSize = 25;

xlabel('\sigma','FontSize', 32)
ylabel('$\sin(\theta_*)$','FontSize', 32,'FontName','Times New Roman','Interpreter','LaTex')

set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , 35              , ...
    'FontName'  , 'Times New Roman');



