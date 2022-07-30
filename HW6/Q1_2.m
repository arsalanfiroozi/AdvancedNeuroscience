%% Define Map
clearvars
clc
close all
addpath(genpath('matlabGiftiCifti'))

% Parameters
num_trials = 100;
learning_rate = 1;
discounting_factor = 1;
beta = 0.5;

M = zeros(15);
M(10,10) = 1;
M(5,5) = -1;
soft_m = zeros(15,15,4); % Down Up Right Left
for i=1:15
    soft_m(1,i,2)= -inf;
end
for i=1:15
    soft_m(15,i,1)= -inf;
end
for i=1:15
    soft_m(i,1,4)= -inf;
end
for i=1:15
    soft_m(i,15,3)= -inf;
end

% figure;
% heatmap(M);
% heatmap(P(:,:,2));

%% Create Demo
% CreateDemo(soft_m,beta,M,discounting_factor,learning_rate,'Demo_beforeTraining.avi');
%% Learning
set(0,'DefaultFigureVisible','off')
X=randi(15,num_trials,1); Y=randi(15,num_trials,1);
for i=1:num_trials
% for i=1:100
    i
    x=X(i); y=Y(i);
    
    path = [];
    decs = [];
    while((x~=5 || y~=5) && (x~=10 || y~=10))
        if(length(path)>100)
            break;
        end
        path = [path; x y];
        r = rand();
        [dec, decision] = mydec(r,[Prob(x,y,1,soft_m,beta) Prob(x,y,2,soft_m,beta) Prob(x,y,3,soft_m,beta) Prob(x,y,4,soft_m,beta)] );
        decs = [decs dec];
        p_d = decision + [x, y];
        
%         % Update Policy
%         soft_m_new = soft_m;
%         for k=1:4
%             t = move(k,[x,y]);
%             xt = t(1);
%             yt = t(2);
%             if(xt==0 || xt==16 || yt==0 || yt==16)
%                 xt=1;
%             elseif k~=dec
% %                 soft_m(x,y,k) = soft_m(x,y,k) - delta * learning_rate * Prob(x,y,k,soft_m,beta) * M(xt,yt);
%                 soft_m_new(x,y,k) = soft_m(x,y,k) - learning_rate * Prob(x,y,k,soft_m,beta) * M(xt,yt);
%             else
%                 soft_m_new(x,y,k) = soft_m(x,y,k) + learning_rate * (1-Prob(x,y,k,soft_m,beta)) * M(xt,yt);
%             end
%         end
%         soft_m = soft_m_new;
%         % Update Value
%         delta = calc_delta(M,[x,y],discounting_factor,soft_m,beta);
%         M(x,y) = M(x,y) + learning_rate*delta;
        
        x = p_d(1);
        y = p_d(2);
    end
%     reward = -(x==5 && y==5) + (x==10 && y==10);
    soft_m_new = soft_m;
    for j=1:size(path,1)
        % Update Value
        delta = calc_delta(M,[path(j,1),path(j,2)],discounting_factor,soft_m,beta);
        M(path(j,1),path(j,2)) = M(path(j,1),path(j,2)) + learning_rate*delta;
        % Update Policy
        for k=1:4
            t = move(k,[path(j,1),path(j,2)]);
            xt = t(1);
            yt = t(2);
            if(xt==0 || xt==16 || yt==0 || yt==16)
                xt=1;
            elseif k~=dec
%                 soft_m(x,y,k) = soft_m(x,y,k) - delta * learning_rate * Prob(x,y,k,soft_m,beta) * M(xt,yt);
                soft_m_new(path(j,1),path(j,2),k) = soft_m(path(j,1),path(j,2),k) - learning_rate * Prob(path(j,1),path(j,2),k,soft_m,beta) * M(xt,yt);
            else
                soft_m_new(path(j,1),path(j,2),k) = soft_m(path(j,1),path(j,2),k) + learning_rate * (1-Prob(path(j,1),path(j,2),k,soft_m,beta)) * M(xt,yt);
            end
        end
    end
    soft_m = soft_m_new;
%     if(mod(i,250)==1)
%         figure;
%         set(gcf,'Color',[1 1 1]);
%         set(gca,'FontName','arial','FontSize',10); % Check this
%         hold on
%         p = [];
%         u = [];
%         for i2=1:15
%             for j2=1:15
%                 t = 0;
%                 for k=1:4
%                     t = t + Prob(i2,j2,k,soft_m,beta) * move(k,[0 0]);
%                 end
%                 t = t ./ norm(t);
%                 p = [p; i2 j2];
%                 u = [u; t];
%             end
%         end
%         quiver(p(:,1),p(:,2),u(:,1),u(:,2));
%         export_fig(['demo/g' num2str(i) '.png'],'-r600');
%     end
end
% writerObj = VideoWriter('Gradient.avi');
% writerObj.FrameRate = 4;
% open(writerObj);
% for j=1:num_trials
%     if(mod(j,250)==1)
%         img =  imread(['demo/g' num2str(j) '.png']);
%     %     frame = im2frame(imresize(img,[393 486]));
%         writeVideo(writerObj, img);
%     end
% end
% close(writerObj);
close all
set(0,'DefaultFigureVisible','on')
heatmap(M);
%% Create Demo
CreateDemo(soft_m,beta,M,discounting_factor,learning_rate,'Demo_afterTraining.avi');
%% Gradients
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
hold on
p = [];
u = [];
for i=1:15
    for j=1:15
        t = 0;
        for k=1:4
            t = t + Prob(i,j,k,soft_m,beta) * move(k,[0 0]);
        end
%         t = t ./ max(abs(t));
        p = [p; i j];
        u = [u; t];
    end
end
quiver(p(:,1),p(:,2),u(:,1),u(:,2),'b');
hold on
scatter(10,10,'filled','r');
scatter(5,5,'filled','black');
title('Transitions');
export_fig('Gradient_Transitions.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
hold on
[fx, fy] = gradient(M);
p = [];
u = [];
for i=1:15
    for j=1:15
        p = [p; i j];
        u = [u; fy(i,j) fx(i,j)];
    end
end
quiver(p(:,1),p(:,2),u(:,1),u(:,2));
hold on
scatter(10,10,'filled','r');
scatter(5,5,'filled','black');
% title('States');
export_fig('Gradient_Values.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
h = heatmap(M);
h.YDisplayData = flipud(h.YDisplayData); 
export_fig('State_Values.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
contour(M,'ShowText','on')
hold on
scatter(10,10,'filled','r');
scatter(5,5,'filled','black');
export_fig('Contour.png','-r600');