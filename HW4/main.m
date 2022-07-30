%% Load Data
clc
clearvars
addpath(genpath('matlabGiftiCifti'))
load('ArrayData/ArrayData.mat');
load('ArrayData/CleanTrials.mat');

Fs = 200;
data_PSD = [];
% Wavelet_data = [];
for i=1:length(chan)
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    for j=1:size(Intersect_Clean_Trials,1)
%         [Wavelet_data(i,j,:,:),f] = cwt(x(:,j)',Fs);
        N = length(x(:,j));
        xdft = fft(x(:,j));
        xdft = xdft(1:floor(N/2)+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/length(x(:,j)):Fs/2;
        data_PSD(i,j,:) = psdx;
    end
end
data_mean_PSD = squeeze(mean(data_PSD,2));
% Wavelet_mean_data = mean(mean(Wavelet_data,2),1);

% Plot PSD
% for i=1:size(data_mean_PSD,1)
%     figure;
%     set(gcf,'Color',[1 1 1]);
%     set(gca,'FontName','arial','FontSize',10); % Check this
%     plot(freq,10*log10(data_mean_PSD(i,:)))
%     grid on
%     title(['PSD of Channel ' num2str(i)])
%     xlabel('Frequency (Hz)')
%     ylabel('Power/Frequency (dB/Hz)')
%     export_fig(['result/PSD_' num2str(i) '.png'],'-r600');
% end
% close all

dominant_freq = ChannelPosition;
Coef = [];
for i=1:size(data_mean_PSD,1)
    i
%     f = figure;
%     set(gcf,'Color',[1 1 1]);
%     set(gca,'FontName','arial','FontSize',10); % Check this
%     plot(freq,log10(data_mean_PSD(i,:)));
%     hold on
%     plot(freq,log10(1./freq));
%     plot(freq,log10(data_mean_PSD(i,:))+log10(freq));
    x = log10(freq);
    y = log10(data_mean_PSD(i,:));
    c = polyfit(x(2:end),y(2:end),1);
    Coef(i,:) = c;
%     plot(freq,polyval(c,log10(freq)),'--');
%     plot(freq,log10(data_mean_PSD(i,:))-polyval(c,log10(freq)));
    y = y - min(y(2:end));
    z = log10(data_mean_PSD(i,:))-polyval(c,log10(freq));
    [ma, mai] = max(z-min(z(2:end))); % Second Approach
%     [ma, mai] = max(log10(data_mean_PSD(i,:))+log10(freq)); % First Approach
    dominant_freq(ChannelPosition == i) = freq(mai);
%     grid on
%     title(['PSD of Channel ' num2str(i)])
%     xlabel('Frequency (Hz)')
%     ylabel('Power/Frequency (dB/Hz)')
%     legend({'PSD' 'Pink Noise, First Approach' 'Cancel Noise by Pink Noise' 'Fitted Line, Second Approach' 'Cancel Noise by the Fitted Line'});
%     export_fig(['result/PSD_PinkCancelled_' num2str(i) '.png'],'-r600');
%     close(f)
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
heatmap(dominant_freq);
title('Dominant Frequencies');
export_fig('DominantFreqs.png','-r600');

save('Dominant_freqs.mat','dominant_freq','Coef');
%% Power Spectogram
clc
clearvars
addpath(genpath('matlabGiftiCifti'))
load('ArrayData/ArrayData.mat');
load('ArrayData/CleanTrials.mat');
load('Dominant_freqs.mat');

Fs = 200;
for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    s_data = [];
    for j=1:size(Intersect_Clean_Trials,1)
%         [s,f,t] = spectrogram(x(:,j),80,60,100,Fs,'yaxis');
        [s,f,t] = stft(x(:,j),Fs);
        s = s(length(s)/2:end,:);
        f = f(length(f)/2:end);
        
        xt = log10(f);
        yt = log10(abs(s));
        for k=1:length(t)
            c = polyfit(xt(2:end),yt(2:end,k),1);
            s(:,k) = 10.^(log10(abs(s(:,k)))-polyval(c,log10(f))); % Approach 2
%             s(:,k) = 10.^(log10(s(:,k))+log10(f)); % Approach 1
        end
        
        s_data(j,:) = reshape(s,[size(s,1)*size(s,2) 1]);
    end
    s_data = mean(s_data,1);
    s_data = reshape(s_data,[size(s,1) size(s,2)]);
    s_data = abs(s_data) .^ 2;
    
    fi = figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    surf(t-1.2, f, 20*log10(abs(s_data)));
    view([0 0 90])
    xlim([t(1)-1.2 t(end)-1.2])
    xlabel('Time (s)')
    ylabel('Freqs. (Hz)')
    title(['Power Spectogram of Channel ' num2str(i)]);
    c = colorbar;
    c.Label.String = 'Power/frequency (dB/Hz)';
    export_fig(['result/PS_C' num2str(i)],'-r600');
    close(fi)
end
close all

Fs = 200;
s_data_d = [];
for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    s_data = [];
    for j=1:size(Intersect_Clean_Trials,1)
%         [s,f,t] = spectrogram(x(:,j),80,60,100,Fs,'yaxis');
        [s,f,t] = stft(x(:,j),Fs);
        s = s(length(s)/2:end,:);
        f = f(length(f)/2:end);
        
        xt = log10(f);
        yt = log10(abs(s));
        for k=1:length(t)
            c = polyfit(xt(2:end),yt(2:end,k),1);
%             s(:,k) = 10.^(log10(abs(s(:,k)))-polyval(c,log10(f))); % Approach 2
            s(:,k) = 10.^(log10(s(:,k))+log10(f)); % Approach 1
        end
        
        s_data(j,:) = reshape(s,[size(s,1)*size(s,2) 1]);
    end
    s_data = mean(s_data,1);
    s_data = reshape(s_data,[size(s,1) size(s,2)]);
    s_data = abs(s_data) .^ 2;
    s_data_d(i,:,:) = s_data;
end
s_data = squeeze(mean(s_data_d,1));
fi = figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    surf(t-1.2, f, 20*log10(abs(s_data)));
    view([0 0 90])
    xlim([t(1)-1.2 t(end)-1.2])
    xlabel('Time (s)')
    ylabel('Freqs. (Hz)')
    title(['Averaged Power Spectogram']);
    c = colorbar;
    c.Label.String = 'Power/frequency (dB/Hz)';
    export_fig('result/Average_PS.png','-r600');
    close(fi)

%% Filter
clc
clearvars
addpath(genpath('matlabGiftiCifti'))
load('ArrayData/ArrayData.mat');
load('ArrayData/CleanTrials.mat');
load('Dominant_freqs.mat');

fd = mean(mean([dominant_freq(1,2:end) dominant_freq(2,2:end) dominant_freq(3,:) dominant_freq(4,:) dominant_freq(5,2:end)]));

Fs = 200;
[b,a] = butter(2,[fd-0.5 fd+0.5]./(Fs/2));
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
freqz(b,a,[],Fs)
title('Filter Specification');
export_fig('Filter.png','-r600');

for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
%     x = mean(x,2);
    x = x(:,10);
    y = filter(b,a,x);
    fi = figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    plot(Time, x);
    hold on
    plot(Time, y);
    legend({'Actual Signal' 'Filtered Signal'})
    xlim([Time(1),Time(end)]);
    xlabel('Time (s)');
    ylabel('LFP Signal (V)');
    export_fig(['result/Filter_C' num2str(i)],'-r600');
    close(fi)
end
close all

%% Phase
clc
clearvars
addpath(genpath('matlabGiftiCifti'))
load('ArrayData/ArrayData.mat');
load('ArrayData/CleanTrials.mat');
load('Dominant_freqs.mat');

fd = mean(mean([dominant_freq(1,2:end) dominant_freq(2,2:end) dominant_freq(3,:) dominant_freq(4,:) dominant_freq(5,2:end)]));

Fs = 200;
[b,a] = butter(2,[fd-0.5 fd+0.5]./(Fs/2));
sigphase_data = [];
for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    sigphase = [];
    for j=1:size(x,2)
        xs = x(:,j);
        y = filter(b,a,xs);
        y = hilbert(y);
        sigphase = [sigphase unwrap(angle(y))];
    end
    sigphase_data(i,:,:) = sigphase;
end

% for i=1:size(sigphase_data,3)
for i=[200]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
        fi = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        heatmap(phase_map)
        caxis([-1, 1]);
        export_fig(['result/S' num2str(j) '.png'],'-r600');
        close(fi);
    end
end
close all

writerObj = VideoWriter('TravellingWave.avi');
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)
    j
    img =  imread(['result/S' num2str(j) '.png']);
    frame = im2frame(img);
    writeVideo(writerObj, frame);
end
close(writerObj);

%% PGD
pgd_data = [];
for i=[50]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end
end

pgd_1 = [];
pgd_2 = [];
pgd_3 = [];
for i=1:641
    x = squeeze(pgd_data(i,:,:,:));
    y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
    pgd_1(i,:,:) = squeeze(y);
end
pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
for i=1:size(pgd_1,1)
    pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
end
pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
for i=1:size(pgd_1,1)
    pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
end
pgd_3 = mean(pgd_3,3);
pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);

PGD = pgd_3 ./ pgd_2';

f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
hist(PGD,20);
title(['PGD Distribution of Trial ' num2str(50)]);
xlabel('PGD');
ylabel('#Count');
export_fig(['PGD_Trial' num2str(50) '.png'],'-r100');
f = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
histogram(PGD(Time>=0),10);
hold on
histogram(PGD(Time<0),10);
title(['PGD Distribution of Trial ' num2str(50)]);
xlabel('PGD');
ylabel('#Count');
legend({'After Cue' 'Before Cue'});
export_fig(['PGD_Trial' num2str(50) '2.png'],'-r100');

%% Direction of Propagation
set(0,'DefaultFigureVisible','off')
for i=[200]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        f = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        quiver(1:10,1:5,GPx,GPy);
        title(['Time: ' num2str(Time(j))]);
        export_fig(['result/GA' num2str(j) '.png'],'-r100');
        close(f)
    end
end
close all

writerObj = VideoWriter('DirectionPropagation_200.avi');
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)
    j
    img =  imread(['result/GA' num2str(j) '.png']);
    frame = im2frame(imresize(img,[393 486]));
    writeVideo(writerObj, frame);
end
close(writerObj);

for i=[1]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        f = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        quiver(1:10,1:5,GPx,GPy);
        title(['Time: ' num2str(Time(j))]);
        export_fig(['result/GA' num2str(j) '.png'],'-r100');
        close(f)
    end
end
close all

writerObj = VideoWriter('DirectionPropagation_1.avi');
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)
    j
    img =  imread(['result/GA' num2str(j) '.png']);
    frame = im2frame(imresize(img,[393 486]));
    writeVideo(writerObj, frame);
end
close(writerObj);
set(0,'DefaultFigureVisible','on')

%% Speed of propagation
speed_data = [];
for i=[100]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        speed_data = [speed_data p];
    end
end
speed_data = speed_data';
speed_data = speed_data(:, [2 5:7 9:48]); 
speed_data_diff = diff(speed_data,1,1) .* Fs;
speed_data_diff = abs(mean(speed_data_diff,2));
speed = speed_data_diff ./ pgd_2(1:end-1)';

fi = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
hist(speed,60);
xlabel('Speed (m/s)');
ylabel('#Count');
xlim([0 1.5]);
box off
export_fig('Speed_Trial_100.png','-r600');
%% Create Demo files, WARNING! ==> This part takes 2 hours to complete
set(0,'DefaultFigureVisible','off')
trial = 100;

pgd_data = [];
for i=[trial]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end
end

pgd_1 = [];
pgd_2 = [];
pgd_3 = [];
for i=1:641
    x = squeeze(pgd_data(i,:,:,:));
    y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
    pgd_1(i,:,:) = squeeze(y);
end
pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
for i=1:size(pgd_1,1)
    pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
end
pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
for i=1:size(pgd_1,1)
    pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
end
pgd_3 = mean(pgd_3,3);
pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);

PGD = pgd_3 ./ pgd_2';

speed_data = [];
for i=[trial]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
%         p = cos(p);
        speed_data = [speed_data p];
    end
end
speed_data = speed_data';
speed_data = speed_data(:, [2 5:7 9:48]); 
speed_data_diff = diff(speed_data,1,1) .* Fs;
speed_data_diff = abs(mean(speed_data_diff,2));
speed = speed_data_diff ./ pgd_2(1:end-1)';

for i=[trial]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)-1
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
        fi = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        heatmap(phase_map)
        caxis([-1, 1]);
        title({['Time: ' num2str(Time(j))] ['Speed: ' num2str(speed(j))] ['PGD: ' num2str(PGD(j))]});
        export_fig(['result/S' num2str(j) '.png'],'-r600');
        close(fi);
    end
end
close all

writerObj = VideoWriter(['Demo_' num2str(trial) '.avi']);
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)-1
    j
    img =  imread(['result/S' num2str(j) '.png']);
    frame = im2frame(img);
    writeVideo(writerObj, frame);
end
close(writerObj);
set(0,'DefaultFigureVisible','on')

set(0,'DefaultFigureVisible','off')
trial = 200;

pgd_data = [];
for i=[trial]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end
end

pgd_1 = [];
pgd_2 = [];
pgd_3 = [];
for i=1:641
    x = squeeze(pgd_data(i,:,:,:));
    y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
    pgd_1(i,:,:) = squeeze(y);
end
pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
for i=1:size(pgd_1,1)
    pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
end
pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
for i=1:size(pgd_1,1)
    pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
end
pgd_3 = mean(pgd_3,3);
pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);

PGD = pgd_3 ./ pgd_2';

speed_data = [];
for i=[trial]
    i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
        j
        phase_map = ChannelPosition;
        p = x(:,j);
%         p = cos(p);
        speed_data = [speed_data p];
    end
end
speed_data = speed_data';
speed_data = speed_data(:, [2 5:7 9:48]); 
speed_data_diff = diff(speed_data,1,1) .* Fs;
speed_data_diff = abs(mean(speed_data_diff,2));
speed = speed_data_diff ./ pgd_2(1:end-1)';

for i=[trial]
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)-1
        j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
        fi = figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10);
        heatmap(phase_map)
        caxis([-1, 1]);
        title({['Time: ' num2str(Time(j))] ['Speed: ' num2str(speed(j))] ['PGD: ' num2str(PGD(j))]});
        export_fig(['result/S' num2str(j) '.png'],'-r600');
        close(fi);
    end
end
close all

writerObj = VideoWriter(['Demo_' num2str(trial) '.avi']);
writerObj.FrameRate = 4;
open(writerObj);
for j=1:size(sigphase_data,2)-1
    j
    img =  imread(['result/S' num2str(j) '.png']);
    frame = im2frame(img);
    writeVideo(writerObj, frame);
end
close(writerObj);
set(0,'DefaultFigureVisible','on')

%% Prefered Direction
clc
trial = 100;

pd_data = [];
for t=1:490
    t
    i = t;
    x = sigphase_data(:,:,i);
    pgd_data = [];
    for j=1:size(sigphase_data,2)
        
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end

    pgd_1 = [];
    pgd_2 = [];
    pgd_3 = [];
    for i=1:641
        x = squeeze(pgd_data(i,:,:,:));
        y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
        pgd_1(i,:,:) = squeeze(y);
    end
    pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
    for i=1:size(pgd_1,1)
        pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
    end
    pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
    for i=1:size(pgd_1,1)
        pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
        pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    end
    pgd_3 = mean(pgd_3,3);
    pgd_3 = squeeze(mean(pgd_3(Time>=0,:),1));
%     pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);
    pd_data(t) = atan(pgd_3(2)/pgd_3(1));
end

fi = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
hist(pd_data,40);
xlabel('Angle');
ylabel('#Count');
title('Prefered Direction');
box off
export_fig('PreferedDirection.png','-r600');
%% Mean Speed
clc
clearvars
addpath(genpath('matlabGiftiCifti'))
load('ArrayData/ArrayData.mat');
load('ArrayData/CleanTrials.mat');
load('Dominant_freqs.mat');

fd = mean(mean([dominant_freq(1,2:end) dominant_freq(2,2:end) dominant_freq(3,:) dominant_freq(4,:) dominant_freq(5,2:end)]));

Fs = 200;
[b,a] = butter(2,[fd-0.5 fd+0.5]./(Fs/2));
sigphase_data = [];
for i=1:length(chan)
    i
    x = chan(i).lfp;
    x = x(:,Intersect_Clean_Trials);
    sigphase = [];
    for j=1:size(x,2)
        xs = x(:,j);
        y = filter(b,a,xs);
        y = hilbert(y);
        sigphase = [sigphase unwrap(angle(y))];
    end
    sigphase_data(i,:,:) = sigphase;
end

speed_data_d = [];
for i=1:490
    speed_data = [];
    t = i
    x = sigphase_data(:,:,i);
    for j=1:size(sigphase_data,2)
%         j
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        speed_data = [speed_data p];
    end
    pgd_data = [];
    for j=1:size(sigphase_data,2)
        
        phase_map = ChannelPosition;
        p = x(:,j);
        p = cos(p);
        for k=1:48
            phase_map(ChannelPosition==k) = p(k);
        end
%         phase_map(1,1) = 0;
%         phase_map(5,1) = 0;
        [GPx, GPy] = gradient(phase_map);
        GPx = GPx  ./ (400*10^(-6));
        GPy = GPy  ./ (400*10^(-6));
        
        pgd_data(j,1,:,:) = GPx;
        pgd_data(j,2,:,:) = GPy;
    end

    pgd_1 = [];
    pgd_2 = [];
    pgd_3 = [];
    for i=1:641
        x = squeeze(pgd_data(i,:,:,:));
        y = sqrt(x(1,:,:).^2+x(2,:,:).^2);
        pgd_1(i,:,:) = squeeze(y);
    end
    pgd_1 = reshape(pgd_1, [size(pgd_1,1) 50]);
    for i=1:size(pgd_1,1)
        pgd_2 = [pgd_2 squeeze(mean(pgd_1(i,~isnan(pgd_1(i,:))),2))];
    end
    pgd_t = reshape(pgd_data, [size(pgd_data,1) 2 50]);
    for i=1:size(pgd_1,1)
        pgd_3(i,1,:) = pgd_t(i,1,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
        pgd_3(i,2,:) = pgd_t(i,2,~isnan(pgd_t(i,1,:)) & ~isnan(pgd_t(i,2,:)));
    end
    pgd_3 = mean(pgd_3,3);
    pgd_3 = squeeze(mean(pgd_3(Time>=0,:),1));
%     pgd_3 = sqrt(pgd_3(:,1).^2+pgd_3(:,2).^2);
    pd_data(t) = atan(pgd_3(2)/pgd_3(1));
    
    speed_data = speed_data';
    speed_data = speed_data(:, [2 5:7 9:48]); 
    speed_data_diff = diff(speed_data,1,1) .* Fs;
    speed_data_diff = abs(mean(speed_data_diff,2));
    speed = speed_data_diff ./ pgd_2(1:end-1)';
%     speed_data_d = [speed_data_d mean(speed)];
    speed_data_d = [speed_data_d speed];
end

speed_data_dr = reshape(speed_data_d, [size(speed_data_d,1)*size(speed_data_d,2) 1]);

fi = figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
hist(speed_data_dr,800);
xlabel('Speed');
ylabel('#Count');
box off
export_fig('Speed1.png','-r600');
xlim([0 1]);
export_fig('Speed2.png','-r600');
