function F = logistic_transition(T, gamma, tau_frac)
% =========================================================================
% DESCRIPTION:
% Logistic transition function
%
%   F(t; gamma, tau) = 1 / (1 + exp(-gamma*(t/T - tau)))
%
% where tau is the midpoint expressed as a fraction of the sample.
% =========================================================================

    t = (1:T)';
    z = -gamma * (t./T - tau_frac);
    F = 1 ./ (1 + exp(z));
end