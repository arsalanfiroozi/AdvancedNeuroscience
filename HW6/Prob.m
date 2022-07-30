function y = Prob(x0,y0,k,m,beta)
    y = 0;
    for i=1:4
        y = y + exp(beta*m(x0,y0,i));
    end
    y = exp(beta*m(x0,y0,k)) / y;
end

