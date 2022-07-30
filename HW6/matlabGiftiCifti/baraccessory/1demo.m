A = {'A' 'B' 'C'};
bars = randi([-400 100],3,5);
bar(categorical(A),bars);

asterisk(2,4,2,3,"Iran",bars);

    xbar = 2;
    ybar = 4;
    xgroup = 2;
    ygroup = 3;
    str_show = "Iran";

    % Values
    nbar = size(bars,2);
    ngroup = size(bars,1);
    
    hold on;
    
    % Pre configuration
    XTicksOld = xticklabels;
    xticksold = xticks;
    height_v = ylim;
    height = height_v(2) - height_v(1);
    xticks(1:ngroup);
    set(gca,'XTickLabel', 1:ngroup);
    
    % Default vars
    errors = zeros(size(bars));
    threshold(1) = height / 125;
    threshold(2) = height * 2 / 125;
    threshold(3) = height * 1.5 / 125;
    threshold(4) = 0;
    
    groupwidth = min(0.8, nbar/(nbar + 1.5));
    %X
    t_x = (1:ngroup) - groupwidth/2 + (2*xbar-1) * groupwidth / (2*nbar);
    x1 = t_x(xgroup);
    t_y = (1:ngroup) - groupwidth/2 + (2*ybar-1) * groupwidth / (2*nbar);
    x2 = t_y(ygroup);
    
    %Y
    y1 = bars(xgroup,xbar);
    y2 = bars(ygroup,ybar);
    
    %P1 & P4
    P1_x = x1;
    P1_y = y1 + (errors(xgroup,xbar) + threshold(1)) * sign(y1);
    P4_x = x2;
    P4_y = y2 + (errors(ygroup,ybar) + threshold(1)) * sign(y2);
    
    %P2 & P3
    P2_x = x1;
    P2_y = P1_y + threshold(2) * sign(P1_y);
    P3_x = x2;
    P3_y = P4_y + threshold(2) * sign(P4_y);
    
    %Correction
    P2_y = max(sign(P2_y)*P2_y,sign(P2_y)*P3_y) * sign(P2_y);
    P3_y = max(sign(P2_y)*P2_y,sign(P2_y)*P3_y) * sign(P2_y);
    
    %Lines
    plot([P1_x P2_x], [P1_y P2_y], 'HandleVisibility','off','Color','black');
    plot([P2_x P3_x], [P2_y P3_y], 'HandleVisibility','off','Color','black');
    plot([P3_x P4_x], [P3_y P4_y], 'HandleVisibility','off','Color','black');
    
    %Text
    text(mean([P2_x, P3_x])+threshold(4), mean([P2_y, P3_y])+threshold(3) * sign(mean([P2_y, P3_y])) - threshold(3) * (sign(mean([P2_y, P3_y])) == -1), str_show, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    
    xticks(xticksold);
    set(gca,'XTickLabel', XTicksOld);
    
    hold off;
