clc
addpath(genpath('matlabGiftiCifti'));
clearvars
close all

B = 0.01;
Sigma = 1;
dt = 0.1;
N = 10000;

decs = [];
RTs = [];
for i=1:N
    i
    [decs(i), RTs(i), x_series] = two_choice_trial(1, -1, Sigma, 0, B, dt);
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
scatter(RTs(decs==1),ones(size(RTs(decs==1))),3,'r','filled');
hold on
scatter(RTs(decs==-1),-1*ones(size(RTs(decs==-1))),3,'b','filled');
legend({'Right' 'Wrong'});
export_fig('Q7.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
RT_sorted = sort(RTs);
d_x = [];
d_y = [];
for i=10:length(RT_sorted)-10
    d_x = [d_x RT_sorted(i)];
    d_y = [d_y sum(decs==-1 & RTs<=RT_sorted(i))/sum(RTs<=RT_sorted(i))*100];
end
plot(d_x,d_y);
% export_fig('Q7_2.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
histogram(RTs(decs==1),300);
hold on
histogram(RTs(decs==-1),300);
title('Inverse Gaussian Distribution');
box off
xlabel('Time (s)');
ylabel('#Number of trials');
legend({'True Decision' 'False Decision'});
% export_fig('Q7_3.png','-r600');
