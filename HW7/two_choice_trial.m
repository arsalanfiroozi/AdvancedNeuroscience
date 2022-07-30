function [d, T, x_series] = two_choice_trial(Tp, Tn, Sigma, x0, B, dt)
    x = x0;
    x_series = [x];
    T = 0;
    while(x<Tp && x>Tn)
        x = x + B * dt + Sigma * normrnd(0,dt);
        x_series = [x_series x];
        T = T + dt;
    end
    d = (x > 0) * 2 - 1;
end

