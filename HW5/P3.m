%% Q2
% Fig 1 Blocking
clearvars

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = 0;
sig = 0.6;
n_Trials = 10;

u = 1;
r = 1;
w_list = w;
sig_list = sig;
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(1) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list sig];
end
r = -1;
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(1) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list sig];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(2,1,1);
plot(0:(2*n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
box off
subplot(2,1,2);
plot(0:(2*n_Trials),sig_list(1,:));
box off
xlabel('#Trials');
ylabel('Sigma');
export_fig('Kalman_Q5.png','-r600');
%% Q6
% Fig 3
clearvars

Gamma = 3;

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(3,1,1);
w0 = [0:-0.02:-1 -3:-0.02:-3.5];
r = w0' + normrnd(0,0.2,length(w0),1);
w0 = w0' + normrnd(0,0.05,length(w0),1);
plot(w0);
hold on
scatter(1:length(w0),r,'x');
ylim([-4 4])
legend({'r(t)' 'w(t)'});
box off
subplot(3,1,2);
plot(w0);
hold on
scatter(1:length(w0),r,'x');
ylim([-4 4])
% export_fig('Kalman_Q5.png','-r600');

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = 0;
sig = 0.6;
u = 1;
w_list = [];
sig_list = [];
beta_list = [];
for i=1:length(r)
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r(i);
    Beta = (y - u*w)^2/(u'*sig*u+V);
    if(Beta > Gamma)
        sig = 1;
    end
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(1) * P;
    w = w_n;
    sig = sig_n;
    beta_list = [beta_list Beta];
    w_list = [w_list w];
    sig_list = [sig_list sig];
end
scatter(1:length(w0),w_list,'o');
legend({'r(t)' 'w(t)' 'w_hat(t)'});
box off
subplot(3,1,3);
plot(beta_list,'-')
hold on
plot([0,length(w0)],[Gamma,Gamma])
legend({'NE' 'Gamma'});
box off
export_fig('Fig3.png','-r600');
%% Q6
% Fig 3
clearvars

Gammas = 0:0.1:15;
mse_list = [];
for iter=1:1000
    iter
    for k=1:length(Gammas)
        Gamma = Gammas(k);
        V = 0.49; % Measure variance
        P = 0.01; % Process variance
        w0 = [0:-0.02:-1 -3:-0.02:-3.5];
        r = w0;
        w0 = w0' + normrnd(0,P,length(w0),1);
        MN = normrnd(0,V,1,length(w0));
        w = 0;
        sig = 0.6;
        u = 1;
        w_list = [];
        sig_list = [];
        beta_list = [];
        for i=1:length(r)
        %     y = r + normrnd(0,V); % This one leads to big difference in result.
            y = r(i) + MN(i);
            Beta = (y - u*w)^2/(u'*sig*u+V);
            if(Beta > Gamma)
                sig = 1;
            end
            w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
            sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(1) * P;
            w = w_n;
            sig = sig_n;
            beta_list = [beta_list Beta];
            w_list = [w_list w];
            sig_list = [sig_list sig];
        end
        mse_list(iter,k) = sum((w_list-r).^2);
    end
end
err = [];
for i=1:size(mse_list,2)
    err = [err std(mse_list(:,i))/sqrt(1000)];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
errorbar(Gammas,mean(mse_list,1),err);
d = mean(mse_list,1);
hold on
plot([Gammas(1) Gammas(end)],ones(1,2)*min(d))
plot([Gammas(1) Gammas(end)],ones(1,2)*d(1))
box off
xlabel('Gamma');
ylabel('MSE');
export_fig('Fig3_d.png','-r600');