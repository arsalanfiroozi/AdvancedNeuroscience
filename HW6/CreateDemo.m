function CreateDemo(soft_m,beta,M,discounting_factor,learning_rate,filename)
    addpath(genpath('matlabGiftiCifti'))
    path = [];
    decs = [];
    decisions = [];
    x=randi(15,1); y=randi(15,1);
    while((x==5 && y==5) || (x==10 && y==10))
        x=randi(15,1); y=randi(15,1);
    end
%     x=15;
%     y=15;
    while((x~=5 || y~=5) && (x~=10 || y~=10))
        if(length(path)>100)
            break;
        end
        path = [path; x y];
        r = rand();
        [dec, decision] = mydec(r,[Prob(x,y,1,soft_m,beta) Prob(x,y,2,soft_m,beta) Prob(x,y,3,soft_m,beta) Prob(x,y,4,soft_m,beta)] );
        decs = [decs dec];
        decisions = [decisions;decision];
        p_d = decision + [x, y];
        x = p_d(1);
        y = p_d(2);
    end
    set(0,'DefaultFigureVisible','off')
    for i=1:length(decisions)-1
        figure;
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontName','arial','FontSize',10); % Check this
        hold on
        scatter(10,10,'filled','r');
        scatter(5,5,'filled','black');
        xlim([1, 15])
        ylim([1, 15])
        axis square
        xticks(1:15);
        yticks(1:15);
        scatter(path(i+1,1),path(i+1,2),'filled','g');
        quiver(path(1:i,1),path(1:i,2),decisions(1:i,1),decisions(1:i,2));
        export_fig(['demo/' num2str(i) '.png'],'-r600');
    end
    close all
    writerObj = VideoWriter(filename);
    writerObj.FrameRate = 4;
    open(writerObj);
    for j=1:length(decisions)-1
        j
        img =  imread(['demo/' num2str(j) '.png']);
        frame = im2frame(imresize(img,[2272 3011]));
        writeVideo(writerObj, frame);
    end
    close(writerObj);
    set(0,'DefaultFigureVisible','on')

end

