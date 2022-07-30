%% Q1
clearvars
clc

addpath(genpath('matlabGiftiCifti'))

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
    
for k=1:3
    load(['S_monkey' num2str(k) '.mat']);
    max_vec = [];
    for i=1:12
        max_vec = [max_vec; sum(S(i).mean_FRs')]
    end

    max_vec = max(max_vec);
    [~, ind] = sort(max_vec);

    data = [];
    for i=1:12
        data = [data max(S(i).mean_FRs(ind(end),:))];
    end
    
    hold on
    plot(0:30:330,data)
    xlabel('Orientation');
    ylabel('Frequency');
end

legend({'Monkey 1' 'Monkey 2' 'Monkey 3'});
export_fig('Q1.png', '-r600');
%% MAP
clearvars
clc
addpath(genpath('matlabGiftiCifti'))

for k=1:3
    load(['data_monkey' num2str(k) '_gratings.mat']);
    indices = 1:150;
    load('removed.mat');
    indices = indices(sum([0;0;removed_neurons{k,1}]==1:length(indices))==0);
    indices = indices(sum([0;0;removed_neurons{k,2}]==1:length(indices))==0);
    
    Map = data.MAP;
    Channels = data.CHANNELS;
    data_mesh = zeros(10,10);
    
    load(['S_monkey' num2str(k) '.mat']);
    
    for i=1:10
        for j=1:10
            [r, c] = find(Channels(:,1)==Map(i,j));
            [c, r] = find(indices == r);
            data_tune = [];
            for q=1:12
                data_tune = [data_tune max(max(S(q).mean_FRs(r,:)))];
            end
            [v, ind] = max(data_tune);
            if(size(ind,1)==0)
                data_mesh(i,j) = nan;
            else
                data_mesh(i,j) = ind;
            end
        end
    end
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    heatmap((data_mesh-1)*30);
    hmap(1:256,1) = linspace(0,1,256);
    hmap(:,[2 3]) = 0.7; %brightness
    huemap = hsv2rgb(hmap);
    colormap(huemap)
    title(['Monkey ' num2str(k)]);
    export_fig(['Q2_' num2str(k) '.png'], '-r600');
end
%% Correlation
clearvars
clc
addpath(genpath('matlabGiftiCifti'))

for k=1:3
    r_sig = [];
    r_noise = [];
    dis = [];
    load(['S_monkey' num2str(k) '.mat']);
    for i=1:size(S(1).mean_FRs,1)
        datai = [];
        for i0=1:12
            datai = [datai sum(S(i0).mean_FRs(i,:))];
        end
        for j=1:i
            dataj = [];
            for i0=1:12
                dataj = [dataj sum(S(i0).mean_FRs(j,:))];
            end
            r = corrcoef(datai,dataj);
            r_sig(i,j) = r(1,2);
            r_sig(j,i) = r(1,2);
        end
    end
    for i=1:size(S(1).mean_FRs,1)
        datai = [];
        for q=1:12
            t = [];
            for t0=1:200
                t = [t; sum(S(q).trial(t0).counts(i,:))];
            end
            datai = [datai; zscore(t)];
        end
        for j=1:i
            dataj = [];
            for q=1:12
                t = [];
                for t0=1:200
                    t = [t; sum(S(q).trial(t0).counts(j,:))];
                end
                dataj = [dataj; zscore(t)];
            end
            r = corrcoef(datai,dataj);
            r_noise(i,j) = r(1,2);
            r_noise(j,i) = r(1,2);
        end
    end
    load(['data_monkey' num2str(k) '_gratings.mat']);
    Map = data.MAP;
    Channels = data.CHANNELS;
    indices = 1:150;
    load('removed.mat');
    indices = indices(sum([0;0;removed_neurons{k,1}]==1:length(indices))==0);
    indices = indices(sum([0;0;removed_neurons{k,2}]==1:length(indices))==0);
    for i=1:size(S(1).mean_FRs,1)
        [xi, yi] = find(Channels(indices(i),1) == Map);
        for j=1:i
            [xj, yj] = find(Channels(indices(j),1) == Map);
            dis(i,j) = sqrt((0.4*xi-0.4*xj)^2+(0.4*yi-0.4*yj)^2);
            dis(j,i) = sqrt((0.4*xi-0.4*xj)^2+(0.4*yi-0.4*yj)^2);
        end
    end
    r_sig = r_sig(tril(ones(size(r_sig)),-1)==1);
    r_noise = r_noise(tril(ones(size(r_noise)),-1)==1);
    dis = dis(tril(ones(size(dis)),-1)==1);
    % First Plot
    G1 = r_sig > 0.5;
    G2 = r_sig > 0 & r_sig <= 0.5;
    G3 = r_sig > -0.5 & r_sig <=0;
    G4 = r_sig <= -0.5;
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    hold on
    
    [Y,E] = discretize(dis(G1),8);
    d = [];
    err = [];
    t = r_noise(G1);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(dis(G2),8);
    d = [];
    err = [];
    t = r_noise(G2);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(dis(G3),8);
    d = [];
    err = [];
    t = r_noise(G3);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(dis(G4),8);
    d = [];
    err = [];
    t = r_noise(G4);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    xlabel('Distance between electrodes (mm)');
    ylabel('Spike count correlation (rsc)');
    legend({'> 0.5' '0 to 0.5' '-0.5 to 0' '< -0.5'});
    title(['Monkey ' num2str(k)]);
    export_fig(['Q3_' num2str(k) '_1.png'], '-r600');
    
    % Second Plot
    G1 = dis <= 1;
    G2 = dis > 1 & dis <=2;
    G3 = dis > 2 & dis <=3;
    G4 = dis > 3 & dis <=4;
    G5 = dis > 5 & dis <=10;
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    hold on
    
    [Y,E] = discretize(r_sig(G1),8);
    d = [];
    err = [];
    t = r_noise(G1);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(r_sig(G2),8);
    d = [];
    err = [];
    t = r_noise(G2);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(r_sig(G3),8);
    d = [];
    err = [];
    t = r_noise(G3);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
    [Y,E] = discretize(r_sig(G4),8);
    d = [];
    err = [];
    t = r_noise(G4);
    for i=1:9
        d = [d mean(t(Y==i))]
        err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
    end
    errorbar(E,d,err)
    
%     [Y,E] = discretize(r_sig(G5),8);
%     d = [];
%     err = [];
%     t = r_noise(G5);
%     for i=1:9
%         d = [d mean(t(Y==i))]
%         err = [err std(t(Y==i))/sqrt(length(t(Y==i)))];
%     end
%     errorbar(E,d,err)
    
    xlabel('Orientation tuning similarity (rsignal)');
    ylabel('Spike count correlation (rsc)');
    legend({'0-1 mm' '1-2 mm' '2-3 mm' '3-4 mm' '5-10 mm'});
    title(['Monkey ' num2str(k)]);
    export_fig(['Q3_' num2str(k) '_2.png'], '-r600');
    % Third Plot
    rs = linspace(-1,1,11);
    ds = linspace(0,4,9);
    data_heat = [];
    for i=2:length(rs)
        for j=2:length(ds)
            sel = (r_sig>rs(i-1) & r_sig<=rs(i)) .* (dis>ds(j-1) & dis<=ds(j));
            data_heat(i-1,j-1) = mean(r_noise(sel==1));
        end
    end
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
%     heatmap(data_heat)
    contourf(data_heat)
%     xlabel('Orientation tuning similarity (rsignal)');
    xticklabels(0:0.5:4.5);
    yticklabels(-1:2/9:1);
    xlabel('Distance between electrodes (mm)');
    ylabel('Spike count correlation (rsc)');
    title(['Monkey ' num2str(k)]);
    colorbar
    export_fig(['Q3_' num2str(k) '_3.png'], '-r600');
end
