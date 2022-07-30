function [rt, p_LIPs, firing_MT, firing_LIP] = lip_activity(MT_p_values,LIP_weights,LIP_threshold, M)
    % Parameters:
    % MT_p_values - a vector with 2 elements, firing probabilities for the
    % excitatory and inhibitory neurons, resp.
    % LIP_weights - a length 2 vector of weighting factors for the evidence
    % from the excitatory (positive) and
    % inhibitory (negative) neurons
    % LIP_threshold - the LIP firing rate that represents the choice threshold criterion
    % use fixed time scale of 1 ms
    dt = 0.001;
    N = [0 0]; % plus is first, minus is second
    rate = 0.0;
    p_LIPs = [];
    LIP_event_times = [];
    t = 0;
    N_LIP = 0;
    firing_MT = [];
    firing_LIP = [];
    while rate < LIP_threshold && t<0.2
        t
        t = t + 0.001;
        dN = [];
        dN(1) = rand() < MT_p_values(1);
        dN(2) = rand() < MT_p_values(2);
        firing_MT = [firing_MT; dN];
        N = N + dN;
        p_LIP = sum(N*LIP_weights);
        p_LIPs = [p_LIPs p_LIP];
        LIP_event = rand() < p_LIP;
        firing_LIP = [firing_LIP LIP_event];
        if(LIP_event)
            LIP_event_times = [LIP_event_times t];
            N_LIP = N_LIP + 1;
        end
        % check LIP mean rate for last M spikes
        rates = 0;
        if(N_LIP > M)
            rate = M/(t - LIP_event_times(N_LIP - M));
        end
        rates = [rates rate];
    end
    rt = t;
end