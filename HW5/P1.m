%% Q1
% Extinction
clearvars
addpath(genpath('matlabGiftiCifti'));

n_Ptrials = 100;
n_Ntrials = 100;
% v => expected reward
% w => weight
% u => presence of stimuli
% d => difference of actual reward and v
% e => learning rate
w_list = [];
e = 0.05;

w = 0;
u = 1;
for i=1:n_Ptrials
    v = w * u;
    w = w + e * (1-v) * u;
    w_list = [w_list w];
end
u = 1;
for i=1:n_Ntrials
    v = w * u;
    w = w + e * (0-v) * u;
    w_list = [w_list w];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:(n_Ntrials+n_Ptrials),w_list);
xlabel('#Trials');
ylabel('W');
box off
hold on
yL = ylim;
plot([n_Ntrials n_Ntrials],[yL(1) yL(2)])
annotation('textarrow',[0.4 0.5],[0.5 0.5],'String','Extinction')
export_fig('Extinction.png','-r600');

%% Q1
% Partial
clearvars
addpath(genpath('matlabGiftiCifti'));

n_Trials = 4000;
% v => expected reward
% w => weight
% u => presence of stimuli
% d => difference of actual reward (r) and v
% e => learning rate
w_list = [];
e = 0.005;

w = 0;
u = 1;
thresholds = 0.25:0.1:1; % Probability of reward
thresholds_Labels = {'0.25' '0.35' '0.45' '0.55' '0.65' '0.75' '0.85' '0.95'}; % Probability of reward
for j=1:length(thresholds)
    threshold = thresholds(j);
    for i=1:n_Trials
        r = rand() < threshold;
        v = w * u;
        w = w + e * (r-v) * u;
        w_list(j,i) = w;
    end
end

w_list_mean = mean(w_list(:,5000:end),2);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:(n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend(thresholds_Labels,'Location', 'eastoutside');
box off
export_fig('Partial.png','-r600');
%     0.2553
%     0.3490
%     0.4450
%     0.5517
%     0.6451
%     0.7606
%     0.8495
%     0.9497
%% Q1
% Blocking
clearvars
addpath(genpath('matlabGiftiCifti'));

n_Trials = 200;
% v => expected reward
% w => weight
% u => presence of stimuli
% d => difference of actual reward (r) and v
% e => learning rate
w1_list = [];
w2_list = [];
e = 0.05;

w1 = 0;
u = 1;
for i=1:n_Trials
    v = w1 * u;
    w1 = w1 + e * (1-v) * u;
    w1_list = [w1_list w1];
end
u = [1 1];
w2 = 0;
w = [w1 w2];
for i=1:n_Trials
    v = [w1 w2] * u';
    w = w + e * (1-v) * u;
    w1_list = [w1_list w(1)];
    w2_list = [w2_list w(2)];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:(2*n_Trials),w1_list);
hold on
plot(n_Trials+1:(2*n_Trials),w2_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
box off
export_fig('Blocking.png','-r600');
%% Q1
% Inhibitory
clearvars
addpath(genpath('matlabGiftiCifti'));

n_Trials = 800;
% v => expected reward
% w => weight
% u => presence of stimuli
% d => difference of actual reward (r) and v
% e => learning rate
w_list = [];
e = 0.05;

w1 = 0;
w2 = 0;
w = [w1;w2];
t = 0.05;
r = 0;
rand_vec = rand(1, n_Trials);
for i=1:n_Trials
    if rand_vec(i)<=0.5
        u = [1;0];
        r = 1;
    else
        u = [1;1];
        r = 0;
    end
    v = u' * w;
    w = w + e * (r-v) * u;
    w_list = [w_list w];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:(n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
box off
export_fig('Inhibitory.png','-r600');
%% Q1
% Overshadow
clearvars
addpath(genpath('matlabGiftiCifti'));

n_Trials = 100;
% v => expected reward
% w => weight
% u => presence of stimuli
% d => difference of actual reward (r) and v
% e1, e2 => learning rate
w_list = [];
e1 = 0.01;
e2 = 0.05;

w1 = 0;
w2 = 0;
w = [w1;w2];
r = 0;
rand_vec = rand(1, n_Trials);
for i=1:n_Trials
    u = [1;1];
    r = 1;
    v = u' * w;
    w = w + [e1; e2] .* (r-v) .* u;
    w_list = [w_list w];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(1:(n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
title(['e1 = ' num2str(e1) ', e2 = ' num2str(e2)]);
box off
export_fig('Overshadow.png','-r600');