function out = cv_supw_test(resid_break, dates_reg_dt, trim_frac, hac_lag, query_k)
% =========================================================================
% DESCRIPTION:
% Computes a one-break SupW test for conditional volatility using the
% Sensier-style dependent variable:
%
%   s_t = sqrt(pi/2) * abs(ehat_t)
%
% where ehat_t are residuals from the conditional-mean break model.
%
% Then tests for a break in the mean of s_t:
%
%   s_t = alpha + delta * I(t > tau) + u_t
%
% using HAC Wald statistics and taking the supremum over trimmed
% candidate break dates.
%
% INPUTS:
%   resid_break  = residuals from CM break model (T x 1)
%   dates_reg_dt = datetime vector aligned with residuals
%   trim_frac    = trimming fraction, e.g. 0.15
%   hac_lag      = HAC bandwidth; if empty, uses default rule
%   query_k      = [optional] local index (1-indexed into s after NaN removal)
%                  at which to also return the Wald stat; used by repartition.
%
% OUTPUT:
%   out              = struct with CV break results
%   out.W_at_query    = Wald stat at query_k (NaN if not requested)
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
    out.supW          = NaN;
    out.break_idx     = NaN;
    out.break_date    = NaT;
    out.mu1           = NaN;
    out.mu2           = NaN;
    out.ratio_post_pre = NaN;
    out.pct_change    = NaN;
    out.trim_n        = NaN;
    out.n_candidates  = NaN;
    out.ssr_nobreak   = NaN;
    out.ssr_best      = NaN;
    out.ssr_gain_pct  = NaN;
    out.cv_proxy      = [];
    out.W_at_query    = NaN;
    out.pval_at_query = NaN;

    if isempty(resid_break) || isempty(dates_reg_dt)
        return;
    end

    ehat = resid_break(:);
    d_use = dates_reg_dt(:);

    keep = ~isnan(ehat);
    ehat = ehat(keep);
    d_use = d_use(keep);

    T = length(ehat);
    out.nobs = T;

    if T < 20
        return;
    end

    % Sensier-style conditional standard deviation proxy
    s = sqrt(pi/2) * abs(ehat);
    out.cv_proxy = s;

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

    % No-break SSR
    mu0 = mean(s);
    ssr0 = sum((s - mu0).^2);
    out.ssr_nobreak = ssr0;

    Wvals = NaN(length(candidates),1);
    mu1_vals = NaN(length(candidates),1);
    mu2_vals = NaN(length(candidates),1);
    ssr_vals = NaN(length(candidates),1);

    for ii = 1:length(candidates)
        k = candidates(ii);

        D = zeros(T,1);
        D((k+1):end) = 1;

        X = [ones(T,1), D];
        beta = X \ s;
        u = s - X * beta;

        ssr_vals(ii) = sum(u.^2);

        XX_inv = inv(X' * X);
        S = zeros(2,2);

        for t = 1:T
            xt = X(t,:)';
            S = S + (u(t)^2) * (xt * xt');
        end

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
        var_delta = Vbeta(2,2);

        if isfinite(var_delta) && var_delta > 0
            Wvals(ii) = (beta(2)^2) / var_delta;
            mu1_vals(ii) = beta(1);
            mu2_vals(ii) = beta(1) + beta(2);
        end
    end

    [supW, idx_best] = max(Wvals);

    % Optional query: return Wald stat at a specific candidate position.
    % Used by repartition to obtain stats/pvals at the final break dates.
    %
    % P-value is 1 - chi2cdf(W, 1):  after repartition the break date is
    % fixed, so W ~ chi2(1) under H0.  SupF table is incorrect here.
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

    out.success = true;
    out.supW = supW;
    out.break_idx = best_k;
    out.break_date = d_use(best_k);
    out.mu1 = mu1_vals(idx_best);
    out.mu2 = mu2_vals(idx_best);
    out.ssr_best = ssr_vals(idx_best);

    if out.mu1 ~= 0
        out.ratio_post_pre = out.mu2 / out.mu1;
        out.pct_change = 100 * (out.mu2 - out.mu1) / out.mu1;
    end

    if ssr0 > 0
        out.ssr_gain_pct = 100 * (ssr0 - out.ssr_best) / ssr0;
    end
end