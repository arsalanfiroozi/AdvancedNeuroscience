clc
addpath(genpath('matlabGiftiCifti'));
clearvars

Ts = 0.5:0.1:10;
B = 0.1;
Sigma = 1;
dt = 0.1;

Err = [];
for j=1:length(Ts)
    j
    T = Ts(j);
    decs = [];
    for i=1:10
        [decs(i), ~] = simple_model(B, Sigma, dt, T);
    end
    Err = [Err sum(decs==-1)/length(decs)*100];
end


figure;
set(gcf,'Color',[1 1 1]);
set(gca,'FontName','arial','FontSize',10); % Check this
plot(Ts,Err)
box off
xlabel('Time');
ylabel('Error %');
export_fig('Q3_10.png','-r600');