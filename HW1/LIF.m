%% a
clearvars
clear

tm = 0.01;
RI = 20/1000;
dt = 0.001;
v = 0;
v_th = 15/1000;

t = 0;
d = [v];
while t<0.1
    v = v + (RI - v)/tm*dt;
    d = [d; v];
    t = t + dt;
    if v > v_th
        v = 0;
    end
end

figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
hold on;
plot(0:dt:0.1,d);
xlabel('time');
export_fig('3_a.png','-r600');
%% c
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));

d = [];
delta_t = 0.00001;
Landa = 100;
t = 0;
while(t<1)
    r = rand(1);
    if r<Landa*delta_t
        d = [d t];
    end
    t = t + delta_t;
end
figure;
subplot(2,1,1);
hold on;
for j=1:size(d,2)
    plot([d(j) d(j)], [0 0.01], 'blue');
end

dt = 0.0001;
t_peak = 0.0015;
% x = 0:dt:1;
% y = 0;
% for i=1:length(d)
%     T = (x-d(i)).*exp(-(x-d(i))/t_peak) .* (x>d(i)) ./ 0.0055;
%     y = y + T;
% end
% plot(x,y)
% % xlim([0 0.2])
% title('I');
spikes_t = 0:dt:1;
spikes = zeros(size(spikes_t));
for i=1:length(d)
    [~,ind] = min(abs(spikes_t-d(i)));
    spikes(ind) = 1;
end
t_kernel = 0:dt:10*t_peak;
kernel_val = kernel(t_kernel,t_peak)/max(kernel(t_kernel,t_peak));

y = 20e-3*conv(spikes, kernel_val, 'same');
plot(spikes_t,y,'red')
title('Input Current');

tm = 0.005;
v = 0;
v_th = 15/1000;
t = 0;
d = [v];
t_spk = [];
while t<1
    if floor(t/dt+1)<=length(y)
        v = v + (y(floor(t/dt+1)) - v)/tm*dt;
    else
        v = v + (0 - v)/tm*dt;
    end
    d = [d; v];
    t = t + dt;
    if v > v_th
        v = 0;
        t_spk = [t_spk; t];
    end
end

subplot(2,1,2);
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot([0 1+dt],[v_th v_th], 'blue');
hold on
plot(0:dt:1+dt,d, 'black');
xlim([0 1])
box off
title('Membrane Voltage');
export_fig('3_c.png','-r600');
%%
clc
clearvars
K = [1:9 10:10:100];
T1 = [0.1:0.1:0.9 1:10]/1000;
CV = [];
for i1=1:length(K)
    i1
    for i2=1:length(T1)
        d = [];
        delta_t = 0.0001;
        Landa = 150;
        t = 0;
        while(t<5)
            r = rand(1);
            if r<Landa*delta_t
                d = [d t];
            end
            t = t + delta_t;
        end
        dt = 0.0001;
        t_peak = 0.0015;
%         x = 0:dt:1;
%         y = 0;
%         for i=1:length(d)
%             T = (x-d(i)).*exp(-(x-d(i))/t_peak) .* (x>d(i)) ./ 0.0055;
%             y = y + T;
%         end
        spikes_t = 0:dt:1;
        spikes = zeros(size(spikes_t));
        for i=1:length(d)
            [~,ind] = min(abs(spikes_t-d(i)));
            spikes(ind) = 1;
        end
        t_kernel = 0:dt:10*t_peak;
        kernel_val = kernel(t_kernel,t_peak)/max(kernel(t_kernel,t_peak));

        y = 20e-3*conv(spikes, kernel_val, 'same');
        tm = T1(i2);
        v = 0;
        v_th = 15/1000;
        t = 0;
        d = [v];
        t_spk = [];
        kf = 0;
        while t<5
            if floor(t/dt+1)<=length(y)
                v = v + (y(floor(t/dt+1)) - v)/tm*dt;
            else
                v = v + (0 - v)/tm*dt;
            end
            d = [d; v];
            t = t + dt;
            if v > v_th
                v = 0;
                kf = kf + 1;
                if(kf==K(i1))
                    t_spk = [t_spk t];
                    kf = 0;
                end
            end
        end
        diff_d = diff(t_spk);
        CV(i1,i2) = std(diff_d)/mean(diff_d);
    end
end
CVD = 0.1:0.1:0.8;
figure;
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
r = 0.01;
for j=1:length(CVD)
    [x,y] = find(CV>CVD(j)-r & CV<CVD(j)+r);
    hold on
    loglog(K(x),T1(y));
end
legend({'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'})
xlabel('K');
ylabel('Tau');
export_fig('3_c2.png','-r600')
%%
clc
clearvars
Mag = [1:9 10:10:100];
TP = [0.1:0.1:0.9 1:10]/1000;
CV = [];
for i1=1:length(Mag)
    i1
    for i2=1:length(TP)
        d = [];
        delta_t = 0.0001;
        Landa = 150;
        t = 0;
        while(t<5)
            r = rand(1);
            if r<Landa*delta_t
                d = [d t];
            end
            t = t + delta_t;
        end
        dt = 0.0001;
        t_peak = 0.0015;
%         x = 0:dt:1;
%         y = 0;
%         for i=1:length(d)
%             T = (x-d(i)).*exp(-(x-d(i))/t_peak) .* (x>d(i)) ./ 0.0055;
%             y = y + T;
%         end
        spikes_t = 0:dt:1;
        spikes = zeros(size(spikes_t));
        for i=1:length(d)
            [~,ind] = min(abs(spikes_t-d(i)));
            spikes(ind) = 1;
        end
        t_kernel = 0:dt:10*TP(i2);
        kernel_val = kernel(t_kernel,t_peak)/max(kernel(t_kernel,t_peak));

        y = Mag(i1)*conv(spikes, kernel_val, 'same');
        tm = 0.005;
        v = 0;
        v_th = 15/1000;
        t = 0;
        d = [v];
        t_spk = [];
        kf = 0;
        while t<5
            if floor(t/dt+1)<=length(y)
                v = v + (y(floor(t/dt+1)) - v)/tm*dt;
            else
                v = v + (0 - v)/tm*dt;
            end
            d = [d; v];
            t = t + dt;
            if v > v_th
                v = 0;
                kf = kf + 1;
                if(kf==5)
                    t_spk = [t_spk t];
                    kf = 0;
                end
            end
        end
        diff_d = diff(t_spk);
        CV(i1,i2) = std(diff_d)/mean(diff_d);
    end
end
CVD = 0.1:0.1:0.8;
figure;
set(gca, 'XScale', 'log', 'YScale', 'log');
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
r = 0.01;
for j=1:length(CVD)
    [x,y] = find(CV>CVD(j)-r & CV<CVD(j)+r);
    hold on
    loglog(Mag(x),TP(y));
end
legend({'0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1.0'})
xlabel('Magnitude');
ylabel('T_peak');
export_fig('4_a.png','-r600');
%% d
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));

d = [];
delta_t = 0.00001;
Landa = 100;
t = 0;
%
dt = 0.0001;
t_peak = 0.0015;
%
tm = 0.005;
v_th = 15/1000;
I = [];
for n=1:100
    t=0;
    d = [];
    while(t<1)
        r = rand(1);
        if r<Landa*delta_t
            d = [d t];
        end
        t = t + delta_t;
    end
    spikes_t = 0:dt:1;
    spikes = zeros(size(spikes_t));
    for i=1:length(d)
        [~,ind] = min(abs(spikes_t-d(i)));
        spikes(ind) = 1;
    end
    t_kernel = 0:dt:10*t_peak;
    kernel_val = kernel(t_kernel,t_peak)/max(kernel(t_kernel,t_peak));
    y = 20e-3*conv(spikes, kernel_val, 'same');
    I(n,:) = y;
end
perc = 0:100;
CV = [];
for n=1:length(perc)
    v = 2*(rand(100,1)>(perc(n)/100))-1;
    Is = v(1) .* I(1,:);
    for i=2:100
        Is = Is + v(i) .* I(i,:);
    end
    y = Is;
    tm = 0.005;
    v = 0;
    v_th = 15/1000;
    t = 0;
    d = [v];
    t_spk = [];
    kf = 0;
    while t<1
        if floor(t/dt+1)<=length(y)
            v = v + (y(floor(t/dt+1)) - v)/tm*dt;
        else
            v = v + (0 - v)/tm*dt;
        end
        d = [d; v];
        t = t + dt;
        if v > v_th
            v = 0;
            kf = kf + 1;
            if(kf==1)
                t_spk = [t_spk t];
                kf = 0;
            end
        end
    end
    diff_d = diff(t_spk);
    CV = [CV std(diff_d)/mean(diff_d)];
end
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(perc,CV);
xlabel('% Percentage of Inhibitory Input Neurons');
ylabel('CV');
box off
export_fig('2_d.png','-r600');
%% e
clear
clearvars
addpath(genpath('matlabGiftiCifti\'));

d = [];
delta_t = 0.00001;
Landa = 100;
t = 0;
%
dt = 0.0001;
t_peak = 0.0015;
%
tm = 0.005;
v_th = 15/1000;
S = [];
for n=1:1000
    t=0;
    d = [];
    while(t<10)
        r = rand(1);
        if r<Landa*delta_t
            d = [d t];
        end
        t = t + delta_t;
    end
    spikes_t = 0:dt:1;
    spikes = zeros(size(spikes_t));
    for i=1:length(d)
        [~,ind] = min(abs(spikes_t-d(i)));
        spikes(ind) = 1;
    end
    S(n,:) = spikes;
end
Pan = (1:100) ./ 1000 /dt;
M = 30;
CV = [];
for i=1:length(Pan)
    i
    d = [];
    for n=Pan(i)+1:size(S,2)
        t = sum(S(:,n-Pan(i):n),2)>0;
        if(sum(t) > M)
            d = [d n];
        end
    end
    timev = diff(d)*dt;
    CV = [CV std(timev)/mean(timev)];
end
figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(Pan.*dt,CV);
xlabel('D');
ylabel('CV');
box off
export_fig('2_e.png','-r600');
