function [t, d] = mydec(r,p)
    t = find(isnan(p));
    d = [1 0];
    if(r<=p(1))
        t = 1;
        d=[1 0];
    elseif(r>p(1) && r<=p(1)+p(2))
        t = 2;
        d=[-1 0];
    elseif(r>p(1)+p(2) && r<=p(1)+p(2)+p(3))
        t = 3;
        d=[0 1];
    elseif(r>p(1)+p(2)+p(3) && r<=p(1)+p(2)+p(3)+p(4))
        t = 4;
        d=[0 -1];
    else
        p
        r
    end
end

