%%
addpath('../export_fig-master');

%% Fig1
figure;
set(gcf,'Color',[1 1 1]);
Data = transpose([45 33 18 37 42 30]);
bar(Data);
set(gca,'XTickLabel', {'A' 'B' 'C' 'D' 'E' 'F'});
sigasterisk(1,1,2,5,"*",Data);
ylim([0 50]);
title("Fig 1");
export_fig('Fig1.png','-r600');

%% Fig2
figure;
set(gcf,'Color',[1 1 1]);
Data = [10 20 30;15 20 5;45 90 10;16 17 13;80 10 30;40 60 10];
bar(Data);
set(gca,'XTickLabel', {'A' 'B' 'C' 'D' 'E' 'F'});
sigasterisk(2,3,3,3,"*",Data);
title("Fig 2");
export_fig('Fig2.png','-r600');

%% Fig3
figure;
set(gcf,'Color',[1 1 1]);
bar(Data);
set(gca,'XTickLabel', {'A' 'B' 'C' 'D' 'E' 'F'});
sigasterisk(1,1,2,3,"**",Data);
sigasterisk(1,3,2,1,"*",Data);
title("Fig 3");
ylim([0 95]);
export_fig('Fig3.png','-r600');

%% Fig4
figure;
set(gcf,'Color',[1 1 1]);
Errors = rand(6,3)*5;
bar(Data);
set(gca,'XTickLabel', {'A' 'B' 'C' 'D' 'E' 'F'});
add_errorbar(Errors, Data);
sigasterisk(2,3,1,1,"*",Data,Errors);
title("Fig 4");
export_fig('Fig4.png','-r600');
