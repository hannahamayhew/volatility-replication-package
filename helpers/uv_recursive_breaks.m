function out = uv_recursive_breaks(v, dates_dt, trim_frac, max_breaks, min_seg_len, min_seg_len_eff_override)
% =========================================================================
% DESCRIPTION:
% Recursively searches for multiple breaks in the mean of a univariate
% series using repeated one-break SupW tests, with:
%
%   - global maximum number of breaks
%   - minimum segment length
%   - 15% trim (or user-specified trim) within each tested segment
%
% Strategy:
%   1. Start with the full sample as one segment
%   2. For each active segment, test for one break
%   3. If significant, record the break and split the segment
%   4. Continue until:
%        - no segment has a significant break
%        - max_breaks total accepted breaks reached
%        - remaining segments are too short
%
% INPUTS:
%   v           = T x 1 series (can contain NaNs)
%   dates_dt    = T x 1 datetime vector
%   trim_frac   = e.g. 0.15
%   max_breaks  = e.g. 5
%   min_seg_len = e.g. 60
%   min_seg_len_eff_override = optional scalar override for the effective
%                              recursive minimum segment length used by
%                              the Bai-Perron h-rule. When empty, the
%                              baseline rule max(ceil(trim_frac*T), min_seg_len)
%                              is used.
%
% OUTPUT:
%   out.break_dates
%   out.break_stats
%   out.break_pvals
%   out.n_breaks
% =========================================================================

    if nargin < 3 || isempty(trim_frac)
        trim_frac = 0.15;
    end
    if nargin < 4 || isempty(max_breaks)
        max_breaks = 5;
    end
    if nargin < 5 || isempty(min_seg_len)
        min_seg_len = 60;
    end
    if nargin < 6
        min_seg_len_eff_override = [];
    end

    out = struct();
    out.break_dates           = NaT(0,1);
    out.break_stats           = [];
    out.break_pvals           = [];
    out.n_breaks              = 0;
    out.repartition_iters     = 0;
    out.repartition_converged = true;
    out.repartition_max_delta = 0;

    % Drop NaNs once up front
    keep = ~isnan(v);
    v = v(keep);
    dates_dt = dates_dt(keep);

    T = length(v);

    % Sensier / Bai-Perron h-rule: minimum segment length is ceil(trim_frac * T).
    % This is sample-relative so that consecutive detected breaks are always
    % at least h = ceil(0.15 * T) observations apart — NOT a fixed month count.
    % The caller-supplied min_seg_len acts only as a hard floor (e.g. for very
    % short series where trim_frac * T < min_seg_len).
    if ~isempty(min_seg_len_eff_override)
        min_seg_len_eff = min_seg_len_eff_override;
    else
        min_seg_len_eff = max(ceil(trim_frac * T), min_seg_len);
    end

    if T < max(min_seg_len_eff, 20)
        return;
    end

    % Segment list: each row is [start_idx, end_idx]
    segments = [1, T];

    break_dates = NaT(0,1);
    break_stats = [];
    break_pvals = [];
    break_idxs  = [];   % global indices into v/dates_dt (needed for repartition)

    while ~isempty(segments) && length(break_dates) < max_breaks

        best_seg_row = NaN;
        best_local_k = NaN;
        best_break_date = NaT;
        best_stat = -Inf;
        best_pval = NaN;

        % Search all current segments for the strongest significant break
        for s = 1:size(segments,1)
            i1 = segments(s,1);
            i2 = segments(s,2);

            seg_len = i2 - i1 + 1;
            if seg_len < max(min_seg_len_eff, 20)
                continue;
            end

            v_seg = v(i1:i2);
            d_seg = dates_dt(i1:i2);

            res = uv_supw_test(v_seg, d_seg, trim_frac, []);

            if ~res.success
                continue;
            end

            pval = hansen_supf_pval_general(res.supW, 1, 0.15);

            if pval < 0.05 && res.supW > best_stat
                best_seg_row = s;
                best_local_k = res.break_idx;
                best_break_date = res.break_date;
                best_stat = res.supW;
                best_pval = pval;
            end
        end

        % Stop if no significant break found in any active segment
        if isnan(best_seg_row)
            break;
        end

        % Record the best currently available break
        break_dates(end+1,1) = best_break_date;
        break_stats(end+1,1) = best_stat;
        break_pvals(end+1,1) = best_pval;

        % Split the winning segment
        i1 = segments(best_seg_row,1);
        i2 = segments(best_seg_row,2);
        k_global = i1 + best_local_k - 1;
        break_idxs(end+1,1) = k_global;    % global index into v — used in repartition

        left_seg  = [i1, k_global];
        right_seg = [k_global+1, i2];

        % Remove the segment we just split
        segments(best_seg_row,:) = [];

        % Add child segments only if long enough to be worth testing later.
        % Uses min_seg_len_eff = ceil(trim_frac * T_full) — the Bai/Bai-Perron h-rule.
        if (left_seg(2) - left_seg(1) + 1) >= max(min_seg_len_eff, 20)
            segments = [segments; left_seg];
        end
        if (right_seg(2) - right_seg(1) + 1) >= max(min_seg_len_eff, 20)
            segments = [segments; right_seg];
        end
    end

    % Sort breaks chronologically before repartition
    [break_dates, idx] = sort(break_dates);
    break_stats = break_stats(idx);
    break_pvals = break_pvals(idx);
    break_idxs  = break_idxs(idx);

    % -----------------------------------------------------------------------
    % REPARTITION: Bai (1997b) iterative one-by-one conditional re-estimation.
    %
    % For each detected break i, hold all other break positions fixed and
    % re-estimate τ_i as the index k ∈ [τ_{i-1}+h, τ_{i+1}-h] minimising:
    %   SSR(v[τ_{i-1}+1 : k]) + SSR(v[k+1 : τ_{i+1}])
    % where h = min_seg_len_eff = ceil(trim_frac * T_full).
    %
    % Cycles over ALL detected breaks repeatedly until the break vector stops
    % changing (integer-grid convergence: max|τ_new − τ_old| = 0) or until
    % MAX_RP_ITER = 100 is reached as a hard stop.
    %
    % After convergence, break_stats and break_pvals are RECOMPUTED at the
    % final repartitioned dates — detection-stage values are NOT retained.
    % For each break i, the Wald stat W is computed at τ_i within segment
    % [τ_{i-1}+1, τ_{i+1}] via uv_supw_test(query_k=τ_i−τ_{i-1}).
    % P-values use 1 - chi2cdf(W, 1): at a fixed known date W ~ chi2(1)
    % under H0 — the SupF table is for suprema over search sets and is
    % incorrect for a fixed-date statistic.
    % -----------------------------------------------------------------------
    rp_iters     = 0;
    rp_converged = true;
    rp_max_delta = 0;

    if ~isempty(break_idxs)
        MAX_RP_ITER = 100;
        h_rp = min_seg_len_eff;
        m_b  = length(break_idxs);
        tau  = sort(break_idxs);

        % Prefix sums for O(1) segment SSR: cs1(j+1) = sum(v(1:j))
        cs1 = [0; cumsum(v)];
        cs2 = [0; cumsum(v .^ 2)];

        rp_converged = false;
        rp_max_delta = Inf;

        while rp_iters < MAX_RP_ITER
            rp_iters = rp_iters + 1;
            tau_prev = tau;

            for ii = 1:m_b
                if ii == 1;   tau_L = 0; else; tau_L = tau(ii-1); end
                if ii == m_b; tau_R = T; else; tau_R = tau(ii+1);  end

                lo = tau_L + h_rp;
                hi = tau_R - h_rp;
                if lo > hi; continue; end

                best_ssr_rp = Inf;
                best_k_rp   = tau(ii);

                for k = lo:hi
                    n1   = k - tau_L;
                    s1   = cs1(k+1)     - cs1(tau_L+1);
                    q1   = cs2(k+1)     - cs2(tau_L+1);
                    ssr1 = max(q1 - s1^2/n1, 0);

                    n2   = tau_R - k;
                    s2   = cs1(tau_R+1) - cs1(k+1);
                    q2   = cs2(tau_R+1) - cs2(k+1);
                    ssr2 = max(q2 - s2^2/n2, 0);

                    if ssr1 + ssr2 < best_ssr_rp
                        best_ssr_rp = ssr1 + ssr2;
                        best_k_rp   = k;
                    end
                end

                tau(ii) = best_k_rp;
            end

            rp_max_delta = max(abs(tau - tau_prev));
            if rp_max_delta == 0
                rp_converged = true;
                break;
            end
        end
        % rp_max_delta = 0 if converged; last-iteration max shift if not.

        % Update break positions from converged tau
        break_idxs  = tau;
        break_dates = dates_dt(break_idxs);

        % Re-sort (repartition preserves order but sort for safety)
        [break_dates, idx2] = sort(break_dates);
        break_idxs  = break_idxs(idx2);
        break_stats = break_stats(idx2);
        break_pvals = break_pvals(idx2);
        tau = break_idxs;

        % -------------------------------------------------------------------
        % RECOMPUTE stats/pvals at final repartitioned dates.
        % For break i: Wald stat at τ_i within segment [τ_{i-1}+1, τ_{i+1}].
        % -------------------------------------------------------------------
        for ii = 1:m_b
            if ii == 1;   tau_L = 0; else; tau_L = tau(ii-1); end
            if ii == m_b; tau_R = T; else; tau_R = tau(ii+1);  end

            seg_v   = v((tau_L+1):tau_R);
            seg_d   = dates_dt((tau_L+1):tau_R);
            k_local = tau(ii) - tau_L;   % 1-indexed position within seg_v

            res_f = uv_supw_test(seg_v, seg_d, trim_frac, [], k_local);
            if ~isnan(res_f.W_at_query)
                break_stats(ii) = res_f.W_at_query;
                break_pvals(ii) = res_f.pval_at_query;
            end
            % On numerical failure (NaN), detection-stage value is retained
            % as a fallback — flagged by the NaN test above.
        end
    end

    out.break_dates           = break_dates;
    out.break_stats           = break_stats;
    out.break_pvals           = break_pvals;
    out.n_breaks              = length(break_dates);
    out.repartition_iters     = rp_iters;
    out.repartition_converged = rp_converged;
    out.repartition_max_delta = rp_max_delta;
end
