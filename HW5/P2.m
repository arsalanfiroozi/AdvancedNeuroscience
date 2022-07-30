%% Q2
% Fig 1 Blocking
clearvars

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = [0; 0];
sig = [0.6 0; 0 0.6];
n_Trials = 10;

u = [1 0];
r = 1;
w_list = w;
sig_list = [sig(1,1);sig(2,2)];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2)]];
end
u = [1 1];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2)]];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(2,1,1);
plot(0:(2*n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
box off
title('Blocking');
hold on
yL = ylim;
plot([n_Trials n_Trials],[yL(1) yL(2)],'--','HandleVisibility','off');
subplot(2,1,2);
plot(0:(2*n_Trials),sig_list(1,:));
hold on
plot(n_Trials:(2*n_Trials),sig_list(2,n_Trials+1:end));
box off
xlabel('#Trials');
ylabel('Sigma');
legend({'S1', 'S2'});
hold on
yL = ylim;
plot([n_Trials n_Trials],[yL(1) yL(2)],'--','HandleVisibility','off');
export_fig('Kalman_Blocking.png','-r600');
%% Q2
% Fig 1 Unblocking
clearvars

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = [0; 0];
sig = [0.6 0; 0 0.6];
n_Trials = 10;

u = [1 0];
r = 1;
w_list = w;
sig_list = [sig(1,1);sig(2,2)];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2)]];
end
u = [1 1];
r = 2;
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2)]];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(2,1,1);
plot(0:(2*n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
box off
title('Unblocking');
hold on
yL = ylim;
plot([n_Trials n_Trials],[yL(1) yL(2)],'--','HandleVisibility','off');
subplot(2,1,2);
plot(0:(2*n_Trials),sig_list(1,:));
hold on
plot(n_Trials:(2*n_Trials),sig_list(2,n_Trials+1:end));
box off
xlabel('#Trials');
ylabel('Sigma');
legend({'S1', 'S2'});
hold on
yL = ylim;
plot([n_Trials n_Trials],[yL(1) yL(2)],'--','HandleVisibility','off');
export_fig('Kalman_Unblocking.png','-r600');
%% Q2
% Fig 1 Drift
clearvars

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = [1; 1];
sig = [0.6 0; 0 0.6];
n_Trials = 20;

u = [1 1];
r = 1;
w_list = w;
sig_list = [sig(1,1);sig(2,2)];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    r = rand() < 0.5;
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2)]];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(0:(n_Trials),w_list);
xlabel('#Trials');
ylabel('W');
legend({'W1', 'W2'});
box off
export_fig('Drift_Unblocking.png','-r600');
%% Q2
% Fig 1 Backwrd Blocking
clearvars

V = 0.49; % Measure variance
P = 0.01; % Process variance
w = [0; 0];
sig = [0.6 0; 0 0.6];
n_Trials = 10;

u = [1 1];
r = 1;
w_list = w;
sig_list = [sig(1,1);sig(2,2);sig(1,2);sig(2,1)];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2);sig(1,2);sig(2,1)]];
end
u = [1 0];
for i=1:n_Trials
%     y = r + normrnd(0,V); % This one leads to big difference in result.
    y = r + V;
    w_n = w + sig * u' ./ (u*sig*u'+V) * (y - u*w);
    sig_n = sig - sig * u' ./ (u*sig*u'+V) * u * sig + eye(2) * P;
    w = w_n;
    sig = sig_n;
    w_list = [w_list w];
    sig_list = [sig_list [sig(1,1);sig(2,2);sig(1,2);sig(2,1)]];
end

% figure;
% set(gcf,'Color',[1 1 1]);
% set(gca,'FontName','arial','FontSize',10); % Check this
% subplot(2,1,1);
% plot(0:(2*n_Trials),w_list);
% xlabel('#Trials');
% ylabel('W');
% legend({'W1', 'W2'});
% box off
% subplot(2,1,2);
% plot(0:(2*n_Trials),sig_list(1,:));
% hold on
% plot(n_Trials:(2*n_Trials),sig_list(2,n_Trials+1:end));
% box off
% xlabel('#Trials');
% ylabel('Sigma');
% legend({'S1', 'S2'});
% export_fig('Kalman_BackwardBlocking.png','-r600');

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
subplot(1,3,1);
w1 = -1:0.02:2;
w2 = -1:0.02:2;
[w1,w2] = meshgrid(w1,w2);
X = [w1(:) w2(:)];
y = mvnpdf(X,w_list(:,1)',[sig_list(1,1) sig_list(3,1); sig_list(4,1) sig_list(2,1)]);
y = reshape(y,length(w2),length(w1));
surf(w1,w2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
axis([-3 3 -3 3 0 0.4])
xlim([-1 2])
ylim([-1 2])
xlabel('w1')
ylabel('w2')
zlabel('Probability Density')
view([0 90]);
shading interp
colormap(gray)
title('t = 1');
hold on
scatter3(w_list(1,1),w_list(2,1),0.3,'*');
axis square
subplot(1,3,2);
w1 = -1:0.02:2;
w2 = -1:0.02:2;
[w1,w2] = meshgrid(w1,w2);
X = [w1(:) w2(:)];
y = mvnpdf(X,w_list(:,9)',[sig_list(1,9) sig_list(3,9); sig_list(4,9) sig_list(2,9)]);
y = reshape(y,length(w2),length(w1));
surf(w1,w2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
axis([-3 3 -3 3 0 0.4])
xlim([-1 2])
ylim([-1 2])
xlabel('w1')
ylabel('w2')
zlabel('Probability Density')
view([0 90]);
shading interp
colormap(gray)
title('t = 9');
hold on
zlim([0 100])
scatter3(w_list(1,9),w_list(2,9),50,'*');
axis square
subplot(1,3,3);
w1 = -1:0.02:2;
w2 = -1:0.02:2;
[w1,w2] = meshgrid(w1,w2);
X = [w1(:) w2(:)];
y = mvnpdf(X,w_list(:,19)',[sig_list(1,19) sig_list(3,19); sig_list(4,19) sig_list(2,19)]);
y = reshape(y,length(w2),length(w1));
surf(w1,w2,y)
caxis([min(y(:))-0.5*range(y(:)),max(y(:))])
axis([-3 3 -3 3 0 0.4])
xlim([-1 2])
ylim([-1 2])
xlabel('w1')
ylabel('w2')
zlabel('Probability Density')
view([0 90]);
shading interp
colormap(gray)
title('t = 19');
hold on
zlim([0 100])
scatter3(w_list(1,19),w_list(2,19),50,'*');
axis square
export_fig('Kalman_Fig2.png','-r600');