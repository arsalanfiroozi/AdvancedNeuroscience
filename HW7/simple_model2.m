function d = simple_model2(B, Sigma, T, x0, dt)
    threshold = normcdf(x0,B*T,sqrt(T*Sigma^2*dt));
    if(rand(1)<threshold)
        d = -1;
    else
        d = +1;
    end
end

