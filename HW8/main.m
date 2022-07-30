%%
clearvars
clc
addpath('DatabaseCode');
showEyeData('./Eye tracking database/Eye tracking database/DATA/hp', './Eye tracking database/Eye tracking database/ALLSTIMULI')

% Modify line 30 of showEyeDataAcrossUsers.m
% datafolder = ['./Eye tracking database/Eye tracking database/DATA/' user];
% showEyeDataAcrossUsers('./Eye tracking database/Eye tracking database/ALLSTIMULI', 6)

%%
clearvars
clc
close all
addpath('DatabaseCode');
addpath(genpath('matlabGiftiCifti'));
addpath(genpath('matlabPyrTools-master'));
addpath(genpath('SaliencyToolbox2.3'));
addpath(genpath('voc-release-3.1-win-master'));
addpath(genpath('LabelMeToolbox-master'));
addpath(genpath('FaceDetect'));
addpath(genpath('gbvs'));
addpath(genpath('JuddSaliencyModel\JuddSaliencyModel'));

FileList = dir('Eye tracking database\Eye tracking database\ALLSTIMULI\');
clear files_social;
% for i = 5:size(FileList,1)
for i = 5:5
    a = FileList(i).name;
    saliencyMap = saliency(['Eye tracking database\Eye tracking database\ALLSTIMULI\' a]);
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    imshow(saliencyMap);
%     export_fig(['Result/Saliency_' a],'-r600');
end

score = rocScoreSaliencyVsFixations(saliencyMap,X,Y,origimgsize)

%%
clearvars
clc
close all
addpath('DatabaseCode');
addpath(genpath('matlabGiftiCifti'));
addpath(genpath('matlabPyrTools-master'));
addpath(genpath('SaliencyToolbox2.3'));
addpath(genpath('voc-release-3.1-win-master'));
addpath(genpath('LabelMeToolbox-master'));
addpath(genpath('FaceDetect'));
addpath(genpath('gbvs'));
addpath(genpath('JuddSaliencyModel\JuddSaliencyModel'));

subs = {'ajs' 'CNG' 'emb' 'ems' 'ff' 'hp' 'jcw' 'jw' 'kae' 'krl' 'po' 'tmj' 'tu' 'ya' 'zb'};

% Data = zeros(15,1003,720,2);
i=4;
m_s = 10000;
% for i=4:size(FileList,1)
FileList = dir(['Eye tracking database\Eye tracking database\DATA\' subs{1}]);
ROCs = zeros(1003,8,2);
%%
clc
% for i=51:size(FileList,1)
for i=51:467
    i
    tic
    d1 = [];
    d2 = [];
    for j=1:15
        data = load(['Eye tracking database\Eye tracking database\DATA\' subs{j} '\' FileList(i).name]);
        name_tmp = FileList(i).name;
        name_tmp = name_tmp(1:end-4);
        data = data.(sprintf("%s", name_tmp));
        data = data.DATA.eyeData;
        if(size(data,1)<720)
            continue;
        end
        d1 = [d1;data(1:(size(data,1)/2),1:2)];
        d2 = [d2;data((size(data,1)/2+1):end,1:2)]; 
        m_s = min(size(data,1),m_s);
    end
    for j=8
        saliencyMap = saliency(['Eye tracking database\Eye tracking database\ALLSTIMULI\' name_tmp '.jpeg'], j);
        sel = d1(:,1) > 0 & d1(:,2)>0;
        ROCs(i-3,j,1) = rocScoreSaliencyVsFixations(saliencyMap, d1(sel==1,1), d1(sel==1,2), size(saliencyMap));
        sel = d2(:,1) > 0 & d2(:,2)>0;
        ROCs(i-3,j,2) = rocScoreSaliencyVsFixations(saliencyMap, d2(sel==1,1), d2(sel==1,2), size(saliencyMap));
    end
    sprintf('TIME:')
    toc
end

%%
clearvars
close all
clc
addpath(genpath('matlabGiftiCifti'));
load('backup.mat')

close all
names = {'SubbandFeatures' 'IttiFeatures' 'ColorFeatures' 'TorralbaSaliency' 'HorizonFeatures' 'ObjectFeatures' 'DistToCenterFeatures'};
for i=1:7
    figure;
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontName','arial','FontSize',10); % Check this
    hist(squeeze(ROCs([1:46 48:464],i,:)),30);
    legend({'First 1.5s' 'Second 1.5s'}, 'Location', 'northwest')
    title({names{i} ['Mean: ' num2str(mean(ROCs([1:49 51:467],i,1))) ', ' num2str(mean(ROCs([1:49 51:467],i,2))) ]})
    export_fig(['Result/' num2str(i) '.png'],'-r600')
end
X = squeeze(ROCs([1:49 51:467],:,:));
err = squeeze(std(X,1)./sqrt(466));
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
errorbar(squeeze(mean(X,1)),err)
legend({'First 1.5s' 'Second 1.5s'}, 'Location', 'northwest')
xticklabels(names);
box off
export_fig(['Result/All.png'],'-r600')
%%
hold on
for i=1:7
    histogram(reshape(ROCs([1:49 51:467],i,:),[2*466 1]),50);
end