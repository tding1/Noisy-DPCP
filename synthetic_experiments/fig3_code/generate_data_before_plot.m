clc
clear

rox = linspace(0,1,400);
rex = linspace(0,0.3,300);

x_val = [];
y_val = [];
t1_val = [];
t2_val = [];
for a = rox
    for b = rex
        if (sqrt(a^2+8)-3*a)^1.5*(sqrt(a^2+8)+a)^0.5 < 32*b
            continue;
        end
        p = [1, 0, a^2-1, 4*a*b, 4*b^2];
        r = roots(p);
        r = sort(r);
        t1 = r(3);
        t2 = r(4);
        
        x_val = [x_val, a];
        y_val = [y_val, b];
        t1_val = [t1_val, t1];
        t2_val = [t2_val, t2];
    end
end

save('x_val.txt', 'x_val', '-ascii')
save('y_val.txt', 'y_val', '-ascii')
save('t1_val.txt', 't1_val', '-ascii')
save('t2_val.txt', 't2_val', '-ascii')

% figure
% c = t1_val;
% scatter(x_val, y_val, 30, c+.001,'o', 'filled')
% % colormap(gray)
% xlabel({'$R_{\mathcal{O}/\widehat{\mathcal{X}}}$'},'fontsize',15, 'Interpreter','latex')
% ylabel({'$R_{\widehat{\mathcal{E}}/\widehat{\mathcal{X}}}$'},'fontsize',15, 'Interpreter','latex')
% alpha(0.75)
% set(gca, ...
%     'LineWidth' , 2                     , ...
%     'FontSize'  , 22              , ...
%     'FontName'  , 'Times New Roman');
% colorbar
% 
% figure
% c = t2_val;
% scatter(x_val, y_val, 30, c+.001,'o', 'filled')
% % colormap(gray)
% xlabel({'$R_{\mathcal{O}/\widehat{\mathcal{X}}}$'},'fontsize',15, 'Interpreter','latex')
% ylabel({'$R_{\widehat{\mathcal{E}}/\widehat{\mathcal{X}}}$'},'fontsize',15, 'Interpreter','latex')
% alpha(0.75)
% set(gca, ...
%     'LineWidth' , 2                     , ...
%     'FontSize'  , 22             , ...
%     'FontName'  , 'Times New Roman');
% colorbar
% 
