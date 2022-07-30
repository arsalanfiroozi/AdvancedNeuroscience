function y = kernel(x,t_peak)
    y = x.*exp(-x/t_peak);
end

