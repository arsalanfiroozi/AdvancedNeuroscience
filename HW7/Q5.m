clc
addpath(genpath('matlabGiftiCifti'));
clearvars

T = 1;
B = 0.1;
Sigma = 1;
dt = 0.1;

decs = [];
for i=1:20
    decs(i) = simple_model2(B, Sigma, T, 0, dt);
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
export_fig(['Q5_B' num2str(B) '.png'],'-r600');
