%% Define Map
clearvars
clc
close all
addpath(genpath('matlabGiftiCifti'))

% Parameters
num_trials = 100;
learning_rates = linspace(0,1,4);
discounting_factors = linspace(0,1,4);
lambda = 0.005;
beta = 0.5;
%% Learning
data = zeros(length(learning_rates));
for k1=1:length(learning_rates)
    k1
    learning_rate = learning_rates(k1);
    for k2=1:length(learning_rates)
        discounting_factor = discounting_factors(k2);
        num_steps = [];
        for q=1:100
            M = zeros(15);
            M(10,10) = 1;
%             M(5,10) = 1;
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
            X=randi(15,num_trials,1); Y=randi(15,num_trials,1);
            for i=1:num_trials
                x=X(i); y=Y(i);

                path = [];
                decs = [];
                while((x~=5 || y~=5) && (x~=10 || y~=10))
%                 while((x~=5 || y~=5) && (x~=5 || y~=10) && (x~=10 || y~=10))
                    if(length(path)>500)
                        break;
                    end
                    path = [path; x y];
                    r = rand();
                    [dec, decision] = mydec(r,[Prob(x,y,1,soft_m,beta) Prob(x,y,2,soft_m,beta) Prob(x,y,3,soft_m,beta) Prob(x,y,4,soft_m,beta)] );
                    decs = [decs dec];
                    p_d = decision + [x, y];

                    x = p_d(1);
                    y = p_d(2);
                end
                soft_m_new = soft_m;
                if(isempty(path))
                    continue;
                end
                delta0 = calc_delta(M,[path(end,1),path(end,2)],discounting_factor,soft_m,beta);
                for j=1:size(path,1)
                    % Update Value
                    delta = calc_delta(M,[path(j,1),path(j,2)],discounting_factor,soft_m,beta) + lambda^(size(path,1)-j)*delta0;
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
                num_steps(q,i) = length(decs);
            end
        end
        num_steps = mean(num_steps(:,end-30:end),1);
        data(k1,k2) = mean(num_steps);
    end
end
%% Gradients
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
h = heatmap(learning_rates, discounting_factors, data);
h.YDisplayData = flipud(h.YDisplayData); 
xlabel('Learning Rate');
ylabel('Discounting Factor');
title(['\lambda: ' num2str(lambda)]);
export_fig(['Q5_1_Lambda' num2str(lambda) '.png'],'-r600');