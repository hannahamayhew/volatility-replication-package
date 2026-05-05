function out = uv_onebreak_search(v, dates_dt, trim_frac)
% =========================================================================
% DESCRIPTION:
% Finds the best single break in the mean of a univariate series v_t using
% a trimmed search over candidate break dates.
%
% Model:
%   v_t = mu1 for t <= tau
%   v_t = mu2 for t >  tau
%
% The chosen break date minimizes the sum of squared residuals (SSR).
%
% -------------------------------------------------------------------------
% INPUTS:
%   v         = T x 1 vector (can contain NaNs)
%   dates_dt  = T x 1 datetime vector aligned with v
%   trim_frac = trimming fraction, e.g. 0.15
%
% OUTPUT:
%   out       = struct with break results
%
% NOTES:
%   - Uses only non-missing observations.
%   - Requires enough usable observations after trimming.
%   - This is a break-location routine, not yet a formal significance test.
% =========================================================================

    % Defaults
    if nargin < 3 || isempty(trim_frac)
        trim_frac = 0.15;
    end

    % Initialize output
    out = struct();
    out.success = false;
    out.nobs = NaN;
    out.trim_frac = trim_frac;
    out.break_idx = NaN;
    out.break_date = NaT;
    out.mu1 = NaN;
    out.mu2 = NaN;
    out.delta = NaN;
    out.ratio_post_pre = NaN;
    out.pct_change = NaN;
    out.ssr_best = NaN;
    out.ssr_nobreak = NaN;
    out.trim_n = NaN;
    out.n_candidates = NaN;

    % Keep only non-missing observations
    keep = ~isnan(v);
    v_use = v(keep);
    d_use = dates_dt(keep);

    T = length(v_use);
    out.nobs = T;

    % Need enough observations to search
    if T < 10
        return;
    end

    % Trimming
    trim_n = ceil(trim_frac * T);
    out.trim_n = trim_n;

    % Candidate break positions:
    % break after observation k, so regime 1 = 1:k, regime 2 = k+1:T
    k_min = trim_n;
    k_max = T - trim_n;

    if k_min >= k_max
        return;
    end

    candidates = k_min:k_max;
    out.n_candidates = length(candidates);

    % No-break SSR
    mu0 = mean(v_use);
    out.ssr_nobreak = sum((v_use - mu0).^2);

    % Search over break dates
    best_ssr = Inf;
    best_k = NaN;
    best_mu1 = NaN;
    best_mu2 = NaN;

    for k = candidates
        mu1 = mean(v_use(1:k));
        mu2 = mean(v_use(k+1:end));

        ssr1 = sum((v_use(1:k) - mu1).^2);
        ssr2 = sum((v_use(k+1:end) - mu2).^2);
        ssr  = ssr1 + ssr2;

        if ssr < best_ssr
            best_ssr = ssr;
            best_k = k;
            best_mu1 = mu1;
            best_mu2 = mu2;
        end
    end

    % Save results
    out.success = true;
    out.break_idx = best_k;
    out.break_date = d_use(best_k);   % last observation in regime 1
    out.mu1 = best_mu1;
    out.mu2 = best_mu2;
    out.delta = best_mu2 - best_mu1;
    out.ssr_best = best_ssr;

    if best_mu1 ~= 0
        out.ratio_post_pre = best_mu2 / best_mu1;
        out.pct_change = 100 * (best_mu2 - best_mu1) / best_mu1;
    end
end