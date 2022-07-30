%%%%
% compile mex files
% mex -v cgf.c nrf/brent.c nrf/frprmn.c nrf/linmin.c nrf/mnbrak.c nrf/nrutil.c -Inrf
%%%%
%% Olshausen
clearvars
clc
close all

addpath(genpath('./'))

VAR_GOAL=0.1;

load('../IMAGES.mat');

for i=1:10
    t = IMAGES(:,:,i);
    IMAGES(:,:,i) = IMAGES(:,:,i)/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
    var(t(:))
end

A = rand(256,192)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));

sparsenet
%% Yale
clearvars
clc
close all

addpath(genpath('./'))

VAR_GOAL=0.1;

filelist1 = dir("..\CroppedYale");
% for i=3:length(filelist1)
for i=3:12
    path = "..\CroppedYale\" + filelist1(i).name;
    filelist = dir(path);
    t = path + "\" + filelist(randi([5 length(filelist)],1,1)).name;
    while(contains(t,"bad"))
        t = path + "\" + filelist(randi([5 length(filelist)],1,1)).name;
    end
    img = imresize(imread(t),[168 168]);
    IMAGES(:,:,i-2) = double(img);
end
close all

for i=1:10
    t = IMAGES(:,:,i);
    t = t - mean(t(:));
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

N=size(IMAGES,1);
M=10;

[fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
rho=sqrt(fx.*fx+fy.*fy);
f_0=0.4*N;
filt=rho.*exp(-(rho/f_0).^4);

for i=1:M
    image=IMAGES(:,:,i);  % you will need to provide get_image
    If=fft2(image);
    imagew=real(ifft2(If.*fftshift(filt)));
    IMAGES1(:,i)=reshape(imagew,N^2,1);
end

IMAGES1=sqrt(0.1)*IMAGES1/sqrt(mean(var(IMAGES1)));
IMAGES = reshape(IMAGES1,N,N,M);

for i=1:10
    t = IMAGES(:,:,i);
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

A = rand(256,100)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));

sparsenet
%% MNIST
clearvars
clc
close all

addpath(genpath('./'))

VAR_GOAL=0.1;

load('../mnist.mat')

IMAGES = test.images;
IMAGES = IMAGES(:,:,1:10);

for i=1:10
    t = IMAGES(:,:,i);
    t = t - mean(t(:));
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

N=size(IMAGES,1);
M=10;

[fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
rho=sqrt(fx.*fx+fy.*fy);
f_0=0.4*N;
filt=rho.*exp(-(rho/f_0).^4);

for i=1:M
    image=IMAGES(:,:,i);  % you will need to provide get_image
    If=fft2(image);
    imagew=real(ifft2(If.*fftshift(filt)));
    IMAGES1(:,i)=reshape(imagew,N^2,1);
end

IMAGES1=sqrt(0.1)*IMAGES1/sqrt(mean(var(IMAGES1)));
IMAGES = reshape(IMAGES1,N,N,M);

for i=1:10
    t = IMAGES(:,:,i);
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

A = rand(256,100)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));

sparsenet
%% 101_ObjectCategories
clearvars
clc
close all

addpath(genpath('./'))

VAR_GOAL=0.1;

filelist1 = dir("..\101_ObjectCategories");
for i=3:length(filelist1)
    path = "..\101_ObjectCategories\" + filelist1(i).name + "\image_0001.jpg";
    img = mean(imread(path),3);
    IMAGES(:,:,i-2) = double(imresize(img,[165 165]));
    imshow(uint8(img));
end
close all

for i=1:10
    t = IMAGES(:,:,i);
    t = t - mean(t(:));
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

N=size(IMAGES,1);
M=10;

[fx, fy]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);
rho=sqrt(fx.*fx+fy.*fy);
f_0=0.4*N;
filt=rho.*exp(-(rho/f_0).^4);

for i=1:M
    image=IMAGES(:,:,i);  % you will need to provide get_image
    If=fft2(image);
    imagew=real(ifft2(If.*fftshift(filt)));
    IMAGES1(:,i)=reshape(imagew,N^2,1);
end

IMAGES1=sqrt(0.1)*IMAGES1/sqrt(mean(var(IMAGES1)));
IMAGES = reshape(IMAGES1,N,N,M);

for i=1:10
    t = IMAGES(:,:,i);
    IMAGES(:,:,i) = t/sqrt(var(t(:))/VAR_GOAL);
    t = IMAGES(:,:,i);
end

A = rand(16,100)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));

sparsenet