function out = select_ar_lag_bic(y, pmax)
% =========================================================================
% DESCRIPTION:
% Select AR lag length p for a univariate series y using BIC on a no-break
% AR(p) model with intercept.
%
% Model:
%   y_t = c + phi_1 y_{t-1} + ... + phi_p y_{t-p} + u_t
%
% INPUTS:
%   y    = T x 1 vector (can contain NaNs)
%   pmax = maximum lag order considered
%
% OUTPUT:
%   out.p_best   = selected lag length
%   out.bic_best = minimum BIC
%   out.T_used   = effective sample size at p_best
% =========================================================================

    out = struct();
    out.p_best = NaN;
    out.bic_best = NaN;
    out.T_used = NaN;

    y = y(:);
    y = y(~isnan(y));

    T = length(y);
    if T < 20
        return;
    end

    bic_vals = NaN(pmax+1,1);
    T_vals   = NaN(pmax+1,1);

    for p = 0:pmax
        if T <= p + 5
            continue;
        end

        Y = y((p+1):end);
        X = ones(length(Y),1);

        for L = 1:p
            X = [X, y((p+1-L):(end-L))];
        end

        beta = X \ Y;
        u = Y - X*beta;

        n = length(Y);
        k = size(X,2);
        sigma2 = mean(u.^2);

        if sigma2 > 0
            bic_vals(p+1) = log(sigma2) + (k * log(n) / n);
            T_vals(p+1) = n;
        end
    end

    [bic_best, idx] = min(bic_vals);

    if isfinite(bic_best)
        out.p_best = idx - 1;
        out.bic_best = bic_best;
        out.T_used = T_vals(idx);
    end
end