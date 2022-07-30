%% Poisson Way 1 - a
clear
clearvars

Landa = 100;
d = [];
t=0;
while(t<1)
    rv = exprnd(1/Landa,1,1);
    t = t + rv;
    d = [d t];
end

figure;
hold on;
for i=1:size(d,2)
    plot([d(i) d(i)], [0 1], 'b');
end
xlim([0 1]);
xlabel('Time (s)');
%% b
clear
clearvars

data = [];
diff_d = [];
for i=1:1000
    Landa = 100;
    d = [];
    t=0;
    while(t<1)
        rv = exprnd(1/Landa,1,1);
        t = t + rv;
        d = [d t];
    end
    data = [data; length(d)];
    diff_d = [diff_d diff(d)];
end

figure;
hist(data,25);
%% c
clear
clearvars

data = [];
diff_d = [];
for i=1:1000
    Landa = 100;
    d = [];
    t=0;
    while(t<1)
        rv = exprnd(1/Landa,1,1);
        t = t + rv;
        d = [d t];
    end
    data = [data; length(d)];
    diff_d = [diff_d diff(d)];
end
% TI = diff(d);
TI = diff_d;
[N, C] = hist(TI,50);
figure;
plot(C,N/1000)
hold on
x_d = 0:0.0001:0.1;
plot(x_d,Landa*exp(-Landa*x_d))
legend({'Simulated' 'Theoretical'});
%% Poisson Way 2
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
for i=1:20
    d = [];
    delta_t = 0.001;
    Landa = 100;
    t = 0;
    while(t<1)
        r = rand(1);
        if r<Landa*delta_t
            d = [d t];
        end
        t = t + delta_t;
    end

    hold on;
    for j=1:size(d,2)
        plot([d(j) d(j)], [i-1 i], 'b');
    end
end
xlabel('Time (s)');
ylabel('Trials');
export_fig('1_a.png','-r600');
%%
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));

data = [];
diff_d = [];
for i=1:10000
    d = [];
    delta_t = 0.0001;
    Landa = 100;
    t = 0;
    while(t<1)
        r = rand(1);
        if r<Landa*delta_t
            d = [d t];
        end
        t = t + delta_t;
    end
    data = [data; length(d)];
    diff_d = [diff_d diff(d)];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
[N, C] = hist(data,200);
bar(C,N/sum(N))
hold on
x_d = 70:130;
for x = x_d
    scatter(x,Landa^x/factorial(x)*exp(-Landa),'r')
end
hold on
% histfit(data,[],'poisson')
title("Distribution of Spikes' Count");
export_fig('1_b.png', '-r600');

TI = diff_d;
[N, C] = hist(TI,100);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(C,N/sum(N))
hold on
x_d = 0:0.0001:0.1;
% plot(x_d,Landa*exp(-Landa*x_d))
legend({'Simulated' 'Theoretical'});
export_fig('1_c.png', '-r600');

sprintf(['CV = ' num2str(std(diff_d)/mean(diff_d))])
%% 2 a
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));
k = 5;

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
for i=1:20
    d = [];
    delta_t = 0.001;
    Landa = 100;
    t = 0;
    kf = 0;
    while(t<1)
        r = rand(1);
        if r<Landa*delta_t
            kf = kf + 1;
            if(kf == k)
                d = [d t];
                kf = 0;
            end
        end
        t = t + delta_t;
    end

    hold on;
    for j=1:size(d,2)
        plot([d(j) d(j)], [i-1 i], 'b');
    end
end
xlabel('Time (s)');
ylabel('Trials');
export_fig('2_a.png','-r600');
%% 2 b,c
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));
k = 5;

data = [];
diff_d = [];
for i=1:10000
    d = [];
    delta_t = 0.0001;
    Landa = 100;
    t = 0;
    kf = 0;
    while(t<1)
        r = rand(1);
        if r<Landa*delta_t
            kf = kf + 1;
            if(kf==1)
                d = [d t];
                kf = 0;
            end
        end
        t = t + delta_t;
    end
    data = [data; length(d)];
    diff_d = [diff_d diff(d)];
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
[N, C] = hist(data,10);
bar(C,N/sum(N))
hold on
% x_d = 70:130;
% for x = x_d
%     scatter(x,Landa^x/factorial(x)*exp(-Landa),'r')
% end
% hold on
% histfit(data,[],'poisson')
title("Distribution of Spikes' Count");
export_fig('2_b.png', '-r600');

TI = diff_d;
[N, C] = hist(TI,100);
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(C,N/sum(N))
hold on
x_d = 0:0.0001:0.1;
% plot(x_d,Landa*exp(-Landa*x_d))
legend({'Simulated' 'Theoretical'});
export_fig('2_c.png', '-r600');

sprintf(['CV = ' num2str(std(diff_d)/mean(diff_d))])
%% With Refractory Period
clear
clearvars
t0 = 0.01;

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10);
for i=1:20
    Landa = 100;
    d = [];
    t=0;
    while(t<1)
        rv = exprnd(1/Landa,1,1);
        t = t + rv + t0;
        d = [d t];
    end
    hold on;
    for j=1:size(d,2)
        plot([d(j) d(j)], [i-1 i], 'b');
    end
    length(d)
end
xlim([0 1]);
xlabel('Time (s)');
ylabel('Trials');
%% g
clearvars
clc

t0 = 1/1000*[1 1 1 0 0 0];
data_p = zeros(3,100);
x_d = linspace(00.0001,30/1000,100);
K = [1 4 51 1 4 51];
for q=1:6
    q
    for k=1:100
        k
        if(x_d(k)>=t0(q))
            data = [];
            for i=1:100
                i
                m = x_d(k)-t0(q);
                d = [];
                t=0;
                kf = 0;
                while(t<5)
%                     t
                    rv = exprnd(m,1,1);
                    kf = kf + 1;
                    t = t + rv + t0(q);
                    if(kf == K(q))
                        d = [d t];
                        kf = 0;
                    end
                end
                data = [data diff(d)];
            end
            data_p(q,k) = std(data) / mean(data);
        end
    end
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(data_p')
legend({'K=1, t0=1 ms' 'K=4, t0=1 ms' 'K=51, t0=1 ms' 'K=1, t0=0 ms' 'K=4, t0=0 ms' 'K=51, t0=0 ms'},'Location','eastoutside');
box off
axis square
xlabel('Delta t (ms)');
ylabel('CV');
ylim([0 1.1]);
export_fig('1_g.png','-r600');