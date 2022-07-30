clc
addpath(genpath('matlabGiftiCifti'));
clearvars
close all

T = 10;
B = 0.1;
Sigma = 1;
dt = 0.1;

N = 100;

traj = [];
decs = [];
for i=1:N
    [decs(i), traj(i,:)] = simple_model(B, Sigma, dt, T);
end

means_th = [];
vars = [];
for i=1:size(traj,2)
    meas(i) = mean(traj(:,i));
    vars(i) = var(traj(:,i));
end
stds = sqrt(vars);

% Theoretical
means_th = 0;
vars_th = 0;
for i=2:size(traj,2)
    means_th(i) = means_th(i-1)+B*dt;
    vars_th(i) = i*Sigma^2*dt^2;
end
stds_th = sqrt(vars_th);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
p = plot(0:dt:10,traj);
for i=1:N
    p(i).Color = [p(i).Color 0.5];
end
hold on
p = plot(0:dt:10,meas);
p.LineWidth = 2;
p.Color = [0 0 0];
p = plot(0:dt:10,meas+stds,'--');
p.LineWidth = 2;
p.Color = [0 0 0];
p = plot(0:dt:10,meas-stds,'--');
p.LineWidth = 2;
p.Color = [0 0 0];

p = plot(0:dt:10,means_th);
p.LineWidth = 2;
p.Color = [1 0 0];
p = plot(0:dt:10,means_th+stds_th,'--');
p.LineWidth = 2;
p.Color = [1 0 0];
p = plot(0:dt:10,means_th-stds_th,'--');
p.LineWidth = 2;
p.Color = [1 0 0];
box off
xlabel('Time');
export_fig('Q4.png','-r600');