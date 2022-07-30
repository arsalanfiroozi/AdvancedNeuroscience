function np = move(r,p)
    if(r==1)
        d=[1 0];
    elseif(r==2)
        d=[-1 0];
    elseif(r==3)
        d=[0 1];
    elseif(r==4)
        d=[0 -1];
    end
    np = p + d;
end

