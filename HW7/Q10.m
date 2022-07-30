clc
clearvars
addpath(genpath('matlabGiftiCifti'))

MT_p_values = [0.15;0.1];
LIP_weights = [0.1;-0.1];
LIP_threshold = 10000;
M = 5;

[rt, p_LIPs, firing_MT, firing_LIP] = lip_activity(MT_p_values,LIP_weights,LIP_threshold,M);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
hold on
for i=1:size(firing_MT,1)
    i
    if(firing_MT(i,1))
        plot([i,i]*0.001,[2,3],'red');
    end
    if(firing_MT(i,2))
        plot([i,i]*0.001,[1,2],'blue');
    end
    if(firing_LIP(i))
        plot([i,i]*0.001,[3,4],'black');
    end
end
yticks([1.5 2.5 3.5]);
yticklabels({'Inhibitory MT Neuron' 'Excitatory MT Neuron' 'LIP Neuron'});
title('Raster Plot');
export_fig('Q10.png','-r600');