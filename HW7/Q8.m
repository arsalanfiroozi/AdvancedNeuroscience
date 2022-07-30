clc
addpath(genpath('matlabGiftiCifti'));
clearvars

B = [0.1 0.3];
Sigma = [1 1];
dt = 0.1;
Tp = [3 2];

[dec, T, x] = race_trial(Tp, Sigma, [0 0], B, dt, Inf);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
% plot(0:dt:T+dt,x);
plot(0:dt:T,x);
hold on
plot([0 T+dt],[Tp(1) Tp(1)],'black');
plot([0 T+dt],[Tp(2) Tp(2)],'black');
legend('Value 1','Value 2','','');
box off
export_fig(['Q8.png'],'-r600');
