%% Parsing Data
clear
clc
clearvars
addpath(genpath('./matlabGiftiCifti'));

load('UnitsData.mat');

conditions = [3 -1;...
              3 +1;...
              6 -1;...
              6 +1;...
              9 -1;...
              9 +1];
Data = [];
N_Trials = [];
for i=1:length(Unit)
    d = Unit(i).Trls;
    Cnd = Unit(i).Cnd(:).Value;
    for j=1:length(Unit(i).Cnd)
        indices = Unit(i).Cnd(j).TrialIdx;
        d_sel = d(indices);
        d_s = [];
        for k=1:length(d_sel)
            d_s = [d_s;d_sel{k}];
        end
        N_Trials(i,j) = length(indices);
        Data{i,j} = d_s;
    end
end

%% Mean PSTH for 6 Conditions
Data_Cond = [];
Count_Cond = sum(N_Trials,1);
x_v = -1.2:0.001:2;
window_len = 0.05;
d_plots = [];
for i=1:size(Data,2)
    d = [];
    for j=1:size(Data,1)
        d = [d;Data{j,i}];
    end
    d = sort(d);
    d_plot = [];
    for j=x_v
        d_plot = [d_plot length(d(d>=(j-window_len) & d<(j)))./ Count_Cond(i) ./ min([window_len,j+1.2])];
    end
    d_plots(i,:) = d_plot;
end


figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(x_v, d_plots);
Ylim = ylim;
hold on;
plot(zeros(1,2),Ylim);
txt = 'Cue \rightarrow';
text(0,max(Ylim)-3,txt,'HorizontalAlignment','right')
plot(ones(1,2)*0.3,Ylim);
txt = 'Delay Period \rightarrow';
text(0.3,max(Ylim)-1,txt,'HorizontalAlignment','right')
plot(ones(1,2)*0.9,Ylim);
txt = 'Target \rightarrow';
text(0.9,max(Ylim)-3,txt,'HorizontalAlignment','right')
xlim([-1.2 2]);
legend({'3 -1', '3 +1', '6 -1', '6 +1', '9 -1', '9 +1'},'Location','southeast');
ylim(Ylim);
xlabel('Time (s)')
ylabel('Freq. (Hz)')
%     title(['PSTH of Condition [' num2str(conditions(i,:)) ']']);
export_fig(['PSTH_Cond_' num2str(conditions(i,1)) num2str(conditions(i,2)) '.png'],'-r600');
%% Some Units PSTH for 6 Conditions
close all
Data_Cond = [];
Sel_Units = [100 200 300 400];
Count_Cond = sum(N_Trials,1);
x_v = -1.2:0.001:2;
window_len = 0.5;
for j=Sel_Units
    d_plots = [];
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        for q=x_v
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.2])];
        end
        d_plots(i,:) = d_plot;
    end
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    plot(x_v, d_plots);
    Ylim = ylim;
    hold on;
    plot(zeros(1,2),Ylim);
    xlim([-1.2 2]);
    legend({'3 -1', '3 +1', '6 -1', '6 +1', '9 -1', '9 +1'},'Location','southeast');
    ylim(Ylim);
    xlabel('Time (s)')
    ylabel('Freq. (Hz)')
    title(['Unit ' num2str(j)]);
    export_fig(['Unit' num2str(j) '_PSTH.png'],'-r600');
end
%% Fit GLM Model for All Units
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
window_len = 0.05;
for j=1:size(Data,1)
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; i];
    end
end
model = fitglm(data_glm,label_glm)
p = coefTest(model)
model = fitglm(data_glm,label_glm(randperm(length(label_glm))));
p = coefTest(model)
%% Fit GLM Model for All Units Location
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
window_len = 0.05;
for j=1:size(Data,1)
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; i];
    end
end
label_glm(label_glm==2) = 1;
label_glm(label_glm==4) = 3;
label_glm(label_glm==6) = 5;
model = fitglm(data_glm,label_glm)
p = coefTest(model)
model = fitglm(data_glm,label_glm(randperm(length(label_glm))));
p = coefTest(model)
%% Fit GLM Model for All Units Location
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
window_len = 0.05;
for j=1:size(Data,1)
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; i];
    end
end
label_glm(label_glm==1) = -1;
label_glm(label_glm==3) = -1;
label_glm(label_glm==5) = -1;
label_glm(label_glm==2) = +1;
label_glm(label_glm==4) = +1;
label_glm(label_glm==6) = +1;
model = fitglm(data_glm,label_glm)
p = coefTest(model)
model = fitglm(data_glm,label_glm(randperm(length(label_glm))));
p = coefTest(model)
%% GLM for Units Seperately Expected Value
clear
clc
clearvars
addpath(genpath('./matlabGiftiCifti'));

load('UnitsData.mat');

conditions = [3 -1;...
              3 +1;...
              6 -1;...
              6 +1;...
              9 -1;...
              9 +1];
x_v = -1.2:0.1:2;
window_len = 0.05;
p_vals = [];
for i=1:length(Unit)
    d = Unit(i).Trls;
    Cnd = Unit(i).Cnd(:).Value;
    data_glm = [];
    label_glm = [];
    for j=1:length(Unit(i).Cnd)
        indices = Unit(i).Cnd(j).TrialIdx;
        d_sel = d(indices);
        for k=1:length(d_sel)
            dd = d_sel{k};
            d_plot = [];
            for q=x_v(2:end)
                d_plot = [d_plot length(dd(dd>=(q-window_len) & dd<(q)))./min([window_len, q+1.200001])];
            end
            data_glm = [data_glm; d_plot];
            label_glm = [label_glm; j];
        end
    end
    label_glm(label_glm==1) = -1;
    label_glm(label_glm==3) = -1;
    label_glm(label_glm==5) = -1;
    label_glm(label_glm==2) = +1;
    label_glm(label_glm==4) = +1;
    label_glm(label_glm==6) = +1;
    model = fitglm(data_glm,label_glm);
    p_vals(i) = coefTest(model);
end
% [p_vals_sorted, idx] = sort(p_vals);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
idx = 1:481;
scatter(idx(p_vals>0.05), p_vals(p_vals>0.05), 5,'b','filled');
hold on
scatter(idx(p_vals<=0.05), p_vals(p_vals<=0.05), 5,'r','filled');
xlabel('Units Index');
ylabel('P-Value');
export_fig('SingleUnits_GLM_PValue.png','-r600');
%% GLM for Units Seperately Location
clear
clc
clearvars
addpath(genpath('./matlabGiftiCifti'));

load('UnitsData.mat');

conditions = [3 -1;...
              3 +1;...
              6 -1;...
              6 +1;...
              9 -1;...
              9 +1];
x_v = -1.2:0.1:2;
window_len = 0.05;
p_vals = [];
for i=1:length(Unit)
    d = Unit(i).Trls;
    Cnd = Unit(i).Cnd(:).Value;
    data_glm = [];
    label_glm = [];
    for j=1:length(Unit(i).Cnd)
        indices = Unit(i).Cnd(j).TrialIdx;
        d_sel = d(indices);
        for k=1:length(d_sel)
            dd = d_sel{k};
            d_plot = [];
            for q=x_v(2:end)
                d_plot = [d_plot length(dd(dd>=(q-window_len) & dd<(q)))./min([window_len, q+1.200001])];
            end
            data_glm = [data_glm; d_plot];
            label_glm = [label_glm; j];
        end
    end
    model = fitglm(data_glm,label_glm);
    p_vals(i) = coefTest(model);
end
% [p_vals_sorted, idx] = sort(p_vals);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
idx = 1:481;
scatter(idx(p_vals>0.05), p_vals(p_vals>0.05), 5,'b','filled');
hold on
scatter(idx(p_vals<=0.05), p_vals(p_vals<=0.05), 5,'r','filled');
xlabel('Units Index');
ylabel('P-Value');
export_fig('SingleUnits_GLM_PValue.png','-r600');
%% GLM for Units Seperately
clear
clc
clearvars
addpath(genpath('./matlabGiftiCifti'));

load('UnitsData.mat');

conditions = [3 -1;...
              3 +1;...
              6 -1;...
              6 +1;...
              9 -1;...
              9 +1];
x_v = -1.2:0.1:2;
window_len = 0.05;
p_vals = [];
for i=1:length(Unit)
    d = Unit(i).Trls;
    Cnd = Unit(i).Cnd(:).Value;
    data_glm = [];
    label_glm = [];
    for j=1:length(Unit(i).Cnd)
        indices = Unit(i).Cnd(j).TrialIdx;
        d_sel = d(indices);
        for k=1:length(d_sel)
            dd = d_sel{k};
            d_plot = [];
            for q=x_v(2:end)
                d_plot = [d_plot length(dd(dd>=(q-window_len) & dd<(q)))./min([window_len, q+1.200001])];
            end
            data_glm = [data_glm; d_plot];
            label_glm = [label_glm; j];
        end
    end
    model = fitglm(data_glm,label_glm);
    p_vals(i) = coefTest(model);
end
% [p_vals_sorted, idx] = sort(p_vals);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
idx = 1:481;
scatter(idx(p_vals>0.05), p_vals(p_vals>0.05), 5,'b','filled');
hold on
scatter(idx(p_vals<=0.05), p_vals(p_vals<=0.05), 5,'r','filled');
xlabel('Units Index');
ylabel('P-Value');
export_fig('SingleUnits_GLM_PValue.png','-r600');
%% PCA Considering All Conditions Altogether
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
data_pca = [];
label_pca = [];
window_len = 0.05;
for j=1:size(Data,1)
%     j
    data_glm = [];
    label_glm = [];
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        t0 = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
            t0 = [t0 i];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; t0];
    end
    data_pca(:,j) = reshape(data_glm,[1 size(data_glm,1)*size(data_glm,2)]);
    label_pca(:,j) = reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)]);
    label_glm0 = label_glm;
%     label_pca = [label_pca reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)])]; 
end

[coeff,score,latent] = pca(data_pca);

Transform = coeff(:,1:3);
Transformed_data = data_pca * Transform;
Transformed_data = reshape(Transformed_data,[size(data_glm,1) size(data_glm,2) 3]);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
for i=1:6
    hold on
    time = 1;
    scatter3(Transformed_data(i,time:end,1),Transformed_data(i,time:end,2),Transformed_data(i,time:end,3),3,'filled');
    view(3);
end
legend({'3 -1', '3 +1', '6 -1', '6 +1', '9 -1', '9 +1'});
xlabel('Dim 1');
ylabel('Dim 2');
zlabel('Dim 3');
title('Dimension Reduction for All Units');
export_fig('PCA_DR.png','-r600');
% label_glm(randperm(length(label_glm)))
%% PCA Considering All Conditions Seperately
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
data_pca = [];
label_pca = [];
window_len = 0.05;
for j=1:size(Data,1)
%     j
    data_glm = [];
    label_glm = [];
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        t0 = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
            t0 = [t0 i];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; t0];
    end
    data_pca(:,j) = reshape(data_glm,[1 size(data_glm,1)*size(data_glm,2)]);
    label_pca(:,j) = reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)]);
    label_glm0 = label_glm;
%     label_pca = [label_pca reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)])]; 
end

labels = label_pca(:,1);
data_transform = zeros([size(data_pca,1),3]);
for i=1:6
    data = data_pca(labels==i,:);
    [coeff,score,latent] = pca(data);
    data_transform(labels==i,:) = data * coeff(:,1:3);
end

Transformed_data = reshape(data_transform,[size(data_glm,1) size(data_glm,2) 3]);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
for i=1:6
    hold on
    time = 250;
    scatter3(Transformed_data(i,time:end,1),Transformed_data(i,time:end,2),Transformed_data(i,time:end,3),3,'filled');
    view(3);
end
legend({'3 -1', '3 +1', '6 -1', '6 +1', '9 -1', '9 +1'});
xlabel('Dim 1');
ylabel('Dim 2');
zlabel('Dim 3');
title('Dimension Reduction for All Units');
export_fig('PCA_DR.png','-r600');
% label_glm(randperm(length(label_glm)))
%% Shuffling
clc
data_glm = [];
label_glm = [];
x_v = -1.2:0.01:2;
data_pca = [];
label_pca = [];
window_len = 0.05;
for j=1:size(Data,1)
%     j
    data_glm = [];
    label_glm = [];
    for i=1:size(Data,2)
        d = Data{j,i};
        d_plot = [];
        t0 = [];
        for q=x_v(2:end)
            d_plot = [d_plot length(d(d>=(q-window_len) & d<(q)))./N_Trials(j,i)./min([window_len, q+1.200001])];
            t0 = [t0 i];
        end
        data_glm = [data_glm; d_plot];
        label_glm = [label_glm; t0];
    end
    data_pca(:,j) = reshape(data_glm,[1 size(data_glm,1)*size(data_glm,2)]);
    label_pca(:,j) = reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)]);
    label_glm0 = label_glm;
%     label_pca = [label_pca reshape(label_glm,[1 size(label_glm,1)*size(label_glm,2)])]; 
end

labels = label_pca(:,1);
labels = labels(randperm(length(labels)));
data_transform = zeros([size(data,1),3]);
for i=1:6
    data = data_pca(labels==i,:);
    [coeff,score,latent] = pca(data);
    data_transform(labels==i,:) = data * coeff(:,1:3);
end

Transformed_data = reshape(data_transform,[size(data_glm,1) size(data_glm,2) 3]);

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
for i=1:6
    hold on
    time = 250;
    scatter3(Transformed_data(i,time:end,1),Transformed_data(i,time:end,2),Transformed_data(i,time:end,3),3,'filled');
    view(3);
end
legend({'3 -1', '3 +1', '6 -1', '6 +1', '9 -1', '9 +1'});
xlabel('Dim 1');
ylabel('Dim 2');
zlabel('Dim 3');
title('Dimension Reduction for All Units');
export_fig('PCA_DR.png','-r600');