clc
clear

file_t = fopen('success.ty','r');
A = fscanf(file_t,'%e');
fclose(file_t);

A = reshape(A,8,8)';

MNratio = 0:0.1:0.7;
noise = 0:0.01:0.07;
imagesc(MNratio,noise,A); colormap(gray(256));  caxis([0 1]);
xlabel('$M/(M+N)$','FontSize',30,'FontName','Times New Roman','Interpreter','LaTex');

set(gca,'YDir','normal')
set(gca, ...
    'LineWidth' , 2                     , ...
    'FontSize'  , 30              , ...
    'FontName'  , 'Times New Roman'     , ...
    'YTick', noise,...
    'XTick', MNratio           );
    xtickformat('%.1f')
    ytickformat('%.2f')
   text(-0.2,0.035,'$\sigma$','FontSize',32,'FontName','Times New Roman','Interpreter','LaTex')

set(gcf, 'Color', 'white');