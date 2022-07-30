clc
clearvars
addpath(genpath('matlabGiftiCifti'))
close all

MT_p_values = [0.1;0.1];
LIP_weights = [0.1 -0.1;-0.1 0.1];
LIP_threshold = [10000 10000];
M = 20;

[rt, p_LIPs, firing_MT, firing_LIP, rates] = lip_activity2(MT_p_values,LIP_weights,LIP_threshold,M);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(2,1,1);
hold on
for i=1:size(firing_MT,1)
    i
    if(firing_MT(i,1))
        plot([i,i]*0.001,[2,3],'red');
    end
    if(firing_MT(i,2))
        plot([i,i]*0.001,[1,2],'blue');
    end
    if(firing_LIP(1,i))
        plot([i,i]*0.001,[3,4],'green');
    end
    if(firing_LIP(2,i))
        plot([i,i]*0.001,[4,5],'black');
    end
end
yticks([1.5 2.5 3.5 4.5]);
yticklabels({'Inhibitory MT Neuron' 'Excitatory MT Neuron' 'LIP Neuron 1' 'LIP Neuron 2'});
title('Raster Plot');
subplot(2,1,2)
plot((1:size(firing_MT,1))*0.001,rates')
box off
xlabel('Time (s)');
ylabel('Firing Rate (Hz)');
legend({'LIP Neuron 1' 'LIP Neuron 2'});
export_fig('Q11.png','-r600');
% 
% figure;
% set(gcf,'Color',[1 1 1]);
% set(gca,'FontName','arial','FontSize',10); % Check this
% plot((1:size(firing_MT,1))*0.001,rates');
% export_fig('Q11.png','-r600');

