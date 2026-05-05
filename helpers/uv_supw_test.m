function out = uv_supw_test(v, dates_dt, trim_frac, hac_lag, query_k)
% =========================================================================
% DESCRIPTION:
% Computes a sup-Wald test for a one-time break in the mean of v_t:
%
%   v_t = alpha + delta * I(t > tau) + u_t
%
% using HAC (Newey-West type) standard errors at each candidate tau, then
% taking the supremum over all trimmed candidate break dates.
%
% -------------------------------------------------------------------------
% INPUTS:
%   v         = T x 1 vector (can contain NaNs)
%   dates_dt  = T x 1 datetime vector aligned with v
%   trim_frac = trimming fraction, e.g. 0.15
%   hac_lag   = HAC bandwidth (e.g. floor(4*(T/100)^(2/9)))
%   query_k   = [optional] local index (1-indexed into v after NaN removal)
%               at which to also return the Wald stat; used by repartition
%               to compute aligned stats at final break dates.
%               If empty or omitted, out.W_at_query = NaN.
%
% OUTPUT:
%   out       = struct with sup-Wald results
%   out.W_at_query    = Wald stat at query_k position (NaN if not requested)
%   out.pval_at_query = 1 - chi2cdf(W, 1)  [fixed-date: chi2(1), not SupF]
% =========================================================================

    if nargin < 3 || isempty(trim_frac)
        trim_frac = 0.15;
    end

    if nargin < 5
        query_k = [];
    end

    out = struct();
    out.success       = false;
    out.nobs          = NaN;
    out.trim_frac     = trim_frac;
    out.hac_lag       = NaN;
    out.W_at_query    = NaN;
    out.pval_at_query = NaN;
    out.supW = NaN;
    out.break_idx = NaN;
    out.break_date = NaT;
    out.alpha = NaN;
    out.delta = NaN;
    out.mu1 = NaN;
    out.mu2 = NaN;
    out.ratio_post_pre = NaN;
    out.pct_change = NaN;
    out.trim_n = NaN;
    out.n_candidates = NaN;

    % Drop missing observations
    keep = ~isnan(v);
    v_use = v(keep);
    d_use = dates_dt(keep);

    T = length(v_use);
    out.nobs = T;

    if T < 20
        return;
    end

    if nargin < 4 || isempty(hac_lag)
        hac_lag = floor(4 * (T / 100)^(2/9));
    end
    out.hac_lag = hac_lag;

    trim_n = ceil(trim_frac * T);
    out.trim_n = trim_n;

    k_min = trim_n;
    k_max = T - trim_n;

    if k_min >= k_max
        return;
    end

    candidates = k_min:k_max;
    out.n_candidates = length(candidates);

    Wvals = NaN(length(candidates),1);
    alpha_vals = NaN(length(candidates),1);
    delta_vals = NaN(length(candidates),1);

    for ii = 1:length(candidates)
        k = candidates(ii);

        D = zeros(T,1);
        D((k+1):end) = 1;

        X = [ones(T,1), D];
        beta = X \ v_use;
        u = v_use - X * beta;

        % HAC variance of beta: (X'X)^(-1) * S * (X'X)^(-1)
        XX_inv = inv(X' * X);
        S = zeros(2,2);

        % lag 0 term
        for t = 1:T
            xt = X(t,:)';
            S = S + (u(t)^2) * (xt * xt');
        end

        % weighted autocovariance terms
        for L = 1:hac_lag
            wL = 1 - L / (hac_lag + 1);
            GammaL = zeros(2,2);

            for t = (L+1):T
                xt = X(t,:)';
                xtL = X(t-L,:)';
                GammaL = GammaL + u(t) * u(t-L) * (xt * xtL');
            end

            S = S + wL * (GammaL + GammaL');
        end

        Vbeta = XX_inv * S * XX_inv;

        % Wald test for H0: delta = 0
        var_delta = Vbeta(2,2);

        if isfinite(var_delta) && var_delta > 0
            Wvals(ii) = (beta(2)^2) / var_delta;
            alpha_vals(ii) = beta(1);
            delta_vals(ii) = beta(2);
        end
    end

    % Sup-W
    [supW, idx_best] = max(Wvals);

    % Optional query: return Wald stat at a specific candidate position.
    % Used by repartition to obtain stats/pvals at the final break dates.
    %
    % P-value is 1 - chi2cdf(W, 1):  after repartition converges the break
    % date is fixed (known), so the Wald stat W is asymptotically chi2(1)
    % under H0 — NOT the supremum distribution.  Using the SupF table here
    % would be incorrect (it is calibrated for the maximum over a search set).
    if ~isempty(query_k) && query_k >= k_min && query_k <= k_max
        qi = query_k - k_min + 1;
        if isfinite(Wvals(qi))
            out.W_at_query    = Wvals(qi);
            out.pval_at_query = 1 - chi2cdf(Wvals(qi), 1);
        end
    end

    if ~isfinite(supW)
        return;
    end

    best_k = candidates(idx_best);
    best_alpha = alpha_vals(idx_best);
    best_delta = delta_vals(idx_best);

    out.success = true;
    out.supW = supW;
    out.break_idx = best_k;
    out.break_date = d_use(best_k);
    out.alpha = best_alpha;
    out.delta = best_delta;
    out.mu1 = best_alpha;
    out.mu2 = best_alpha + best_delta;

    if out.mu1 ~= 0
        out.ratio_post_pre = out.mu2 / out.mu1;
        out.pct_change = 100 * (out.mu2 - out.mu1) / out.mu1;
    end
end