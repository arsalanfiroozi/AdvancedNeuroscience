clc
addpath(genpath('matlabGiftiCifti'));
clearvars

T = 1;
B = 1;
Sigma = 1;
dt = 0.1;

decs = [];
for i=1:20
    [decs(i), ~] = simple_model(B, Sigma, dt, T);
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
bar([-0.5 0.5],[sum(decs==-1) sum(decs==1)])
title(['Decision Distribution for B =' num2str(B)]);
xticks([-0.5 0.5]);
xlim([-1 1]);
xticklabels({'No Go' 'Go'});
box off
export_fig(['Q1_B' num2str(B) '.png'],'-r600');

%% Part 2

Bs = [1 10 0 0.1 -1];
T = 10;
decs = [];
x_s = [];
for i=1:5
    B = Bs(i);
    [decs(i), x_s(i,:)] = simple_model(B, Sigma, dt, T);
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(0:dt:T,x_s')
xlim([0 T]);
% xticklabels({'No Go' 'Go'});
box off
legend(arrayfun(@num2str, Bs, 'UniformOutput', 0),'Location','northwest');
xlabel('Time');
export_fig(['Q2.png'],'-r600');