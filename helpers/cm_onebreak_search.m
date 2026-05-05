function out = cm_onebreak_search(y, dates_dt, p, trim_frac)
% =========================================================================
% DESCRIPTION:
% Estimates a one-break AR(p) conditional-mean model by searching over
% break dates and choosing the break that minimizes SSR.
%
% Model:
%   y_t = c1 + phi11 y_{t-1} + ... + phi1p y_{t-p}          for t <= tau
%   y_t = c2 + phi21 y_{t-1} + ... + phi2p y_{t-p} + u_t    for t >  tau
%
% Equivalent regression form:
%   y_t = c + sum phi_i y_{t-i}
%       + d0*I(t>tau) + sum d_i y_{t-i}*I(t>tau) + u_t
%
% INPUTS:
%   y         = T x 1 vector (can contain NaNs)
%   dates_dt  = T x 1 datetime vector aligned with y
%   p         = fixed AR lag order
%   trim_frac = trimming fraction, e.g. 0.15
%
% OUTPUT:
%   out       = struct with best break results
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
    out.break_idx_reg = NaN;
    out.break_date = NaT;
    out.ssr_best = NaN;
    out.ssr_nobreak = NaN;
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

    trim_n = ceil(trim_frac * Treg);
    out.trim_n = trim_n;

    k_min = trim_n;
    k_max = Treg - trim_n;

    if k_min >= k_max
        return;
    end

    candidates = k_min:k_max;
    out.n_candidates = length(candidates);

    % No-break model
    beta0 = Xbase \ Y;
    u0 = Y - Xbase * beta0;
    ssr0 = sum(u0.^2);

    best_ssr = Inf;
    best_k = NaN;
    best_beta = [];

    for k = candidates
        D = zeros(Treg,1);
        D((k+1):end) = 1;

        X = [Xbase, D, Xbase(:,2:end).*D];
        beta = X \ Y;
        u = Y - X*beta;
        ssr = sum(u.^2);

        if ssr < best_ssr
            best_ssr = ssr;
            best_k = k;
            best_beta = beta;
            best_u = u;
        end
    end

    out.success = true;
    out.break_idx_reg = best_k;
    out.break_date = d_reg(best_k);
    out.ssr_best = best_ssr;
    out.ssr_nobreak = ssr0;
    out.beta_nobreak = beta0;
    out.beta_break = best_beta;
    out.resid_nobreak = u0;
    out.resid_break = best_u;

    if ssr0 > 0
        out.ssr_gain_pct = 100 * (ssr0 - best_ssr) / ssr0;
    end
end