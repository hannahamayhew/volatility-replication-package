function out = cm_supw_test(y, dates_dt, p, trim_frac, hac_lag)
% =========================================================================
% DESCRIPTION:
% Computes a sup-Wald test for a one-time break in the conditional mean of
% an AR(p) model:
%
%   y_t = c + phi_1 y_{t-1} + ... + phi_p y_{t-p}
%       + d0*I(t>tau) + d1*y_{t-1}*I(t>tau) + ... + dp*y_{t-p}*I(t>tau) + u_t
%
% using HAC (Newey-West style) covariance estimation at each candidate
% break date, and taking the supremum of the Wald statistic over all
% trimmed candidate dates.
%
% INPUTS:
%   y         = T x 1 vector (can contain NaNs)
%   dates_dt  = T x 1 datetime vector aligned with y
%   p         = fixed AR lag order
%   trim_frac = trimming fraction, e.g. 0.15
%   hac_lag   = HAC bandwidth; if empty, uses default rule
%
% OUTPUT:
%   out       = struct with sup-Wald results
% =========================================================================

    if nargin < 4 || isempty(trim_frac)
        trim_frac = 0.15;
    end

    out = struct();
    out.success = false;
    out.nobs_raw = NaN;
    out.nobs_reg = NaN;
    out.p = p;
    out.trim_frac = trim_frac;
    out.hac_lag = NaN;
    out.m_restr = NaN;
    out.supW = NaN;
    out.break_idx_reg = NaN;
    out.break_date = NaT;
    out.ssr_nobreak = NaN;
    out.ssr_best = NaN;
    out.ssr_gain_pct = NaN;
    out.trim_n = NaN;
    out.n_candidates = NaN;
    out.beta_nobreak = [];
    out.beta_break = [];
    out.resid_nobreak = [];
    out.resid_break = [];

    keep = ~isnan(y);
    y_use = y(keep);
    d_use = dates_dt(keep);

    out.nobs_raw = length(y_use);

    if length(y_use) < p + 20
        return;
    end

    % Build AR regression sample
    Y = y_use((p+1):end);
    Xbase = ones(length(Y),1);
    for L = 1:p
        Xbase = [Xbase, y_use((p+1-L):(end-L))];
    end
    d_reg = d_use((p+1):end);

    Treg = length(Y);
    out.nobs_reg = Treg;

    if nargin < 5 || isempty(hac_lag)
        hac_lag = floor(4 * (Treg / 100)^(2/9));
    end
    out.hac_lag = hac_lag;

    trim_n = ceil(trim_frac * Treg);
    out.trim_n = trim_n;

    k_min = trim_n;
    k_max = Treg - trim_n;

    if k_min >= k_max
        return;
    end

    candidates = k_min:k_max;
    out.n_candidates = length(candidates);

    % Number of restrictions under H0: no break in conditional mean
    m = p + 1;
    out.m_restr = m;

    % No-break model
    beta0 = Xbase \ Y;
    u0 = Y - Xbase * beta0;
    ssr0 = sum(u0.^2);

    Wvals = NaN(length(candidates),1);
    beta_break_store = cell(length(candidates),1);
    resid_break_store = cell(length(candidates),1);
    ssr_store = NaN(length(candidates),1);

    for ii = 1:length(candidates)
        k = candidates(ii);

        D = zeros(Treg,1);
        D((k+1):end) = 1;

        X = [Xbase, D, Xbase(:,2:end).*D];
        beta = X \ Y;
        u = Y - X*beta;
        ssr = sum(u.^2);

        % HAC covariance
        XX_inv = inv(X' * X);

        % lag 0 (vectorized: sum_t u(t)^2 * X(t,:)'*X(t,:))
        Xu = X .* u;
        S  = Xu' * Xu;

        % higher lags (vectorized: sum_t u(t)*u(t-L) * X(t,:)'*X(t-L,:)')
        for L = 1:hac_lag
            wL    = 1 - L / (hac_lag + 1);
            Xu_t  = X((L+1):end,:) .* u((L+1):end);
            Xu_tL = X(1:end-L,:)   .* u(1:end-L);
            GammaL = Xu_t' * Xu_tL;
            S = S + wL * (GammaL + GammaL');
        end

        Vbeta = XX_inv * S * XX_inv;

        % Restriction matrix for break parameters: [d0 d1 ... dp] = 0
        q = size(Xbase,2); % = p+1
        R = [zeros(m,q), eye(m)];
        Rb = R * beta;
        RVRT = R * Vbeta * R';

        if rcond(RVRT) > 1e-12
            Wvals(ii) = Rb' * (RVRT \ Rb);
        end

        beta_break_store{ii} = beta;
        resid_break_store{ii} = u;
        ssr_store(ii) = ssr;
    end

    [supW, idx_best] = max(Wvals);

    if ~isfinite(supW)
        return;
    end

    best_k = candidates(idx_best);

    out.success = true;
    out.supW = supW;
    out.break_idx_reg = best_k;
    out.break_date = d_reg(best_k);
    out.ssr_nobreak = ssr0;
    out.ssr_best = ssr_store(idx_best);
    out.beta_nobreak = beta0;
    out.beta_break = beta_break_store{idx_best};
    out.resid_nobreak = u0;
    out.resid_break = resid_break_store{idx_best};

    if ssr0 > 0
        out.ssr_gain_pct = 100 * (ssr0 - out.ssr_best) / ssr0;
    end
end