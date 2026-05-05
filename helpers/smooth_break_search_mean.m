function out = smooth_break_search_mean(y, dates_dt, trim_frac, gamma_grid)
% =========================================================================
% DESCRIPTION:
% Fits a smooth-transition mean-shift model
%
%   y_t = a + b * F(t; gamma, tau) + u_t
%
% where F is the logistic transition function.
%
% Searches over:
%   - tau in the trimmed region
%   - gamma over a supplied grid
%
% and chooses the combination minimizing SSR.
%
% INPUTS:
%   y          = T x 1 series (may contain NaNs)
%   dates_dt   = T x 1 datetime vector aligned with y
%   trim_frac  = trimming fraction, e.g. 0.15
%   gamma_grid = vector of gamma values, e.g. [1 2 5 10 20 50 100]
%
% OUTPUT:
%   out        = struct with best smooth-transition fit
% =========================================================================

    if nargin < 3 || isempty(trim_frac)
        trim_frac = 0.15;
    end
    if nargin < 4 || isempty(gamma_grid)
        gamma_grid = [1 2 5 10 20 50 100];
    end

    out = struct();
    out.success = false;
    out.nobs = NaN;
    out.gamma_best = NaN;
    out.tau_frac_best = NaN;
    out.midpoint_date = NaT;
    out.ssr_best = NaN;
    out.ssr_nobreak = NaN;
    out.ssr_gain_pct = NaN;
    out.a = NaN;
    out.b = NaN;
    out.pre_mean = NaN;
    out.post_mean = NaN;
    out.ratio_post_pre = NaN;
    out.pct_change = NaN;

    keep = ~isnan(y);
    y_use = y(keep);
    d_use = dates_dt(keep);

    T = length(y_use);
    out.nobs = T;

    if T < 20
        return;
    end

    trim_n = ceil(trim_frac * T);
    tau_idx_grid = trim_n:(T-trim_n);

    if isempty(tau_idx_grid)
        return;
    end

    mu0 = mean(y_use);
    ssr0 = sum((y_use - mu0).^2);
    out.ssr_nobreak = ssr0;

    best_ssr = Inf;

    for g = 1:length(gamma_grid)
        gamma = gamma_grid(g);

        for kk = 1:length(tau_idx_grid)
            k = tau_idx_grid(kk);
            tau_frac = k / T;

            F = logistic_transition(T, gamma, tau_frac);
            X = [ones(T,1), F];

            beta = X \ y_use;
            u = y_use - X*beta;
            ssr = sum(u.^2);

            if ssr < best_ssr
                best_ssr = ssr;
                best_gamma = gamma;
                best_tau_frac = tau_frac;
                best_beta = beta;
                best_k = k;
            end
        end
    end

    out.success = true;
    out.gamma_best = best_gamma;
    out.tau_frac_best = best_tau_frac;
    out.midpoint_date = d_use(best_k);
    out.ssr_best = best_ssr;
    out.a = best_beta(1);
    out.b = best_beta(2);
    out.pre_mean = out.a;
    out.post_mean = out.a + out.b;

    if out.pre_mean ~= 0
        out.ratio_post_pre = out.post_mean / out.pre_mean;
        out.pct_change = 100 * (out.post_mean - out.pre_mean) / out.pre_mean;
    end

    if ssr0 > 0
        out.ssr_gain_pct = 100 * (ssr0 - best_ssr) / ssr0;
    end
end