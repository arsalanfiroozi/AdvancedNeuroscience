function [d, T, x_series] = race_trial(Tp, Sigma, x0, B, dt, Tm)
    x = x0;
    x_series = [x];
    T = 0;
    while(x(1)<Tp(1) && x(2)<Tp(2) && T<Tm)
        x = x + B .* dt + Sigma .* [normrnd(0,dt) normrnd(0,dt)];
        x_series = [x_series;x];
        T = T + dt;
    end
    d = (x(1) > x(2)) * 2 - 1;
end

