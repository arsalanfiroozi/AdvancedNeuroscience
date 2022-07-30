%% Q3
clearvars
clc
close all

VAR_GOAL=0.1;

vidObj = VideoReader('../BIRD.avi');
for i=1:10
    vidFrame = imresize(readFrame(vidObj),[288 288]);
    imshow(vidFrame)
    pause(1/vidObj.FrameRate);
    IMAGES(:,:,i) = mean(vidFrame,3);
end

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
%%
A = rand(256,100)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));
%%
sparsenet

save('Q3.mat');
%%
clearvars
close all
load('Q3.mat');
addpath(genpath('../matlabGiftiCifti'))
display_network(A,S_var);

set(0,'DefaultFigureVisible','off')
k = 1;
Data = [];
pathes_num = [3*18+4 3*18+10 6*18+5 6*18+17 13*18+5 13*18+15];
while hasFrame(vidObj)
    vidFrame = mean(imresize(readFrame(vidObj),[288 288]),3);
    
    t = vidFrame;
    t = t - mean(t(:));
    vidFrame = t/sqrt(var(t(:))/VAR_GOAL);
    
    image=vidFrame;  % you will need to provide get_image
    If=fft2(image);
    imagew=real(ifft2(If.*fftshift(filt)));
    IMAGES1=reshape(imagew,N^2,1);
    IMAGES1=sqrt(0.1)*IMAGES1/sqrt(mean(var(IMAGES1)));
    vidFrame = reshape(IMAGES1,N,N);
    
    t = vidFrame;
    vidFrame = t/sqrt(var(t(:))/VAR_GOAL);
    
    batch_size = (size(IMAGES,2)/sqrt(size(A,1)))^2;
    X = [];
    for i=1:sqrt(batch_size)
        for j=1:sqrt(batch_size)
            r = (i-1)*sqrt(size(A,1))+1;
            c = (j-1)*sqrt(size(A,2))+1;
            X(:,(i-1)*sqrt(batch_size)+j) = reshape(vidFrame(r:r+sz-1,c:c+sz-1),L,1);
        end
    end
    S=cgf_fitS(A,X,noise_var,beta,sigma,tol);
    [x y] = max(S);
    y = reshape(y, [1 1]*sqrt(length(y)));
    h = heatmap(y);
    h.YDisplayData = flipud(h.YDisplayData); 
    export_fig(['../Result/' num2str(k) '.png'],'-r600');
    close all
    for i=pathes_num
        figure;
        plot(S(:,i))
        ylim([-1.5 1.5])
        export_fig(['../Result/' num2str(i) '_' num2str(k) '.png'],'-r600');   
    end
    Data(k,:,:) = S;
    k = k + 1;
    close all
end
set(0,'DefaultFigureVisible','on')

for i=pathes_num
    writerObj = VideoWriter(['../Coef_Changes_' num2str(i) '.avi']);
    writerObj.FrameRate = 4;
    open(writerObj);
    for j=1:108
        img =  imread(['../Result/' num2str(i) '_' num2str(j) '.png']);
        %     frame = im2frame(imresize(img,[393 486]));
        writeVideo(writerObj, img);
    end
    close(writerObj);
end
writerObj = VideoWriter(['../bird_biggest_basis.avi']);
writerObj.FrameRate = 4;
open(writerObj);
for j=1:108
    img =  imread(['../Result/' num2str(j) '.png']);
        frame = im2frame(imresize(img,[2272 2986]));
    writeVideo(writerObj, frame);
end
close(writerObj);
%%
close all
figure;
hold on
for i=1:108
    i
    for j=1:324
        view(3);
        scatter3(i*ones(size(squeeze(Data(i,:,j)))),j*ones(size(squeeze(Data(i,:,j)))),Data(i,:,j))
    end
end
%%
close all
figure;
hold on
[X,Y] = meshgrid(1:324,1:108);
Z = [];
for i=1:size(X,1)
    for j=1:size(X,2)
        [~, Z(i,j)] = max(Data(Y(i,j),:,X(i,j)));
    end
end
surf(X,Y,Z)