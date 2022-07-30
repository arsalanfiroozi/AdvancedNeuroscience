function delta = calc_delta(M,pos,discounting_factor,soft_m,beta)
    x = pos(1);
    y = pos(2);
%     delta = M(x,y)/sum(sum(abs(M)));
    delta = 0;
    t = 0;
    if x~=15
        t = t + discounting_factor * Prob(x,y,1,soft_m,beta) * M(x+1,y);
    end
    if x~=1
        t = t + discounting_factor * Prob(x,y,2,soft_m,beta) * M(x-1,y);
    end
    if y~=15
        t = t + discounting_factor * Prob(x,y,3,soft_m,beta) * M(x,y+1);
    end
    if y~=1
        t = t + discounting_factor * Prob(x,y,4,soft_m,beta) * M(x,y-1);
    end
    delta = delta + t - M(x,y);
end
