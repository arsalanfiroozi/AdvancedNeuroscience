function [d, x_series] = simple_model(B, Sigma, dt, T)
    x = 0;
    x_series = [x];
    for i=0:dt:T-dt
        x = x + B * dt + Sigma * normrnd(0,dt);
        x_series = [x_series x];
    end
    d = (x > 0) * 2 - 1;
end

